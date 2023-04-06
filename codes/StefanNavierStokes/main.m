% SIMPLE algorithm for the incompressible Navier-Stokes equations
%  - advection is explicit (MUSCL scheme) 
%  - viscous terms are implicit  (BTCS scheme)
%  - heat condubtion is implicit (BTCS scheme)
% We solve the momentum and the internal energy equation
% The Boussinesq approximation is used for buyoncy forces

clear;
close all;
% clc;

global Nx Ny dx dy dt nu  uLid g beta T0 lambda uWall vWall TBC xc yc;
global Tlake Tair;
% physical parameters
nu = 1e-2;      % kinematic viscosity
uLid =  1;       % velocity of the lid (top boundary)
uWall = 0;      % u velocity at the walls
vWall = 0;      % v velocity at the walls
g = 9.81;       % accelaration due to gravity
T0 = 0;         % reference temperature
TBC = 2;        % Boundary condition temperature
lambda = 1e-3;  % thermal diffusivity
beta = 0;       % thermal expansion coefficient = drho/dT
Tlake = 4;
Tair  = 4;

% computational parameters
time = 0;
tend = 4;

xL = 0;
xR = 1;
yL = 0;
yR = 1;
Nx = 64;
Ny = 64;
dx = (xR - xL)/Nx;
dy = (yR - yL)/Ny;

xc = linspace(xL+dx/2,xR-dx/2,Nx);    % x coordinates of the cell-centers
yc = linspace(yL+dy/2,yR-dy/2,Ny);    % y coordinates of the cell-centers
x = linspace(xL,xR,Nx+1);             % x coordinates of the cell-edges
y = linspace(yL,yR,Ny+1);             % y coordinates of the cell-edges
CFL = 0.9;                              % CFL number < 1
dt_accel = 1e-3;                        % a number between [0,1], 1 corresponds to the FTCS time step, 0 to fully BTCS

% initial  conditions
u = zeros(Nx+1,Ny);                 % x component of the velocity vector (blue in the lectures) +1 because cell faces
v = zeros(Nx,Ny+1);                 % y component of the velocity vector (green in the lectures)
p = zeros(Nx,Ny);                   % pressure -> cell centre
T = zeros(Nx,Ny);                   % temperature -> cell centre
Tc= zeros(Nx+1,Ny+1);               % temperature -> cell corners
rhs =zeros(Nx,Ny);                  % right hand-side for the pressure Poisson equation

for i = 1:Nx
    for j = 1:Ny
       T(i,j) = T0+4-0*exp(-0.5*( (xc(i) - 0.8 )^2 + (yc(j) - 0.75)^2)/0.1^2 ); %cold bubble
%        T(i,j) = T0+exp(-0.5*( (xc(i) - 0.2 )^2 + (yc(j) - 0.25)^2)/0.1^2 ); %hot bubble
    end
end

% Read data for the reference solution from files:
% data are taken from Sverdrup et al. "Highly parallelisable simulations 
% of time-dependent viscoplastic fluid flow with strubtured adaptive mesh 
% refinement", Phys Fluids 2018;30(9):093102. doi: 10.1063/1.5049202 . 
% http://arxiv.org/ abs/1803.00417 .
uref = readmatrix('ldc_Re100_u_vs_y.dat');
vref = readmatrix('ldc_Re100_v_vs_x.dat');

fig1 = figure;
set(gcf,'Position',[2549 743 1280 495])

% TIME LOOP
nmax = 1000000;
for n=1:nmax

    umax = max(max(abs(u)));
    vmax = max(max(abs(v)));

    % impose the stability restriction on the time step (advection + diffusion)
    % restriction by convective (advection,umax,vamx) and parabolic part (diffusion,nu)
    if max(umax,vmax) < CFL*dx*100/tend
        dt = tend/100;
    else
        dt = CFL/( umax/dx + vmax/dy + dt_accel*2*nu*(1/dx^2 + 1/dy^2) );
    end% 
%     if (dt > 0.1)
%         dt = 0.1;
%     end
    if( time + dt > tend)
        dt = tend - time;
    end
    if (time >= tend)
        break;
    end
    
    % STEP #1
    pstar = p;  % initial guess for the pressure

    % STEP #2 Explicit:
    % compute the nonlinear convective and diffusive terms
    [rhsustar,rhsvstar] = MomentumConvectionDiffusion(u,v,T);
    % add the grad(pstar) to the velocity
    [rhsustar,rhsvstar] = MomentumPressureGrad(rhsustar,rhsvstar,pstar);

    % STEP #3 Implicit:
    Tc = bary2corners(T,Tlake,Tlake,Tlake,Tair);
    [ustar,erru,ku] = CG_u_vel(rhsustar,T,Tc);
    [vstar,errv,kv] = CG_v_vel(rhsvstar,T,Tc);

    % compute the right hand-side for the Poisson equation
    rhs = (ustar(2:Nx+1,:) - ustar(1:Nx,:))/(dt*dx) + ...
          (vstar(:,2:Ny+1) - vstar(:,1:Ny))/(dt*dy);
    [pprime,errp,kp] = ConjGradOpt_p(rhs);  % use the matrix free Conjugate Gradient method to solve the pressure Poisson equation
 
    % STEP #4 Div-free velocity correction
    [u,v] = DivFreeVelocityCorrection(u,v,ustar,vstar,pprime);
    divuv = (u(2:Nx+1,:) - u(1:Nx,:))/dx + (v(:,2:Ny+1) - v(:,1:Ny))/dy;

    % update pressure adding the correction
    p = pstar + pprime;

    % time update
    time = time + dt;

    % velocity at the cell barycenters (only for the vizualization!)
    ub = 0.5*( u(2:Nx+1,:) + u(1:Nx,:) );
    vb = 0.5*( v(:,2:Ny+1) + v(:,1:Ny) );
    
    subplot(1,3,1)
    hold off
    s = surf(xc,yc,sqrt(ub.^2+vb.^2)','EdgeColor','none','FaceColor','interp');
%     s = surf(xc(2:Nx-1),yc(2:Ny-1),abs(divuv(2:Nx-1,2:Ny-1))','EdgeColor','none','FaceColor','interp');
    view([0 90]) % position the camera, vision angle
    hold on 
    quiver(xc,yc,ub',vb','k');
    title(strcat('Current time = ',num2str(time)))
    xlabel('x [m]')
    ylabel('y [m]')
    %Z0 = get(s,'ZData');
    %set(s,'ZData',Z0 - 10)
    axis square;
    colorbar;


    % plot 1D cuts:
    subplot(1,3,2)
    plot(ub(Nx/2,:)',yc,'LineWidth',1.5)
    hold on
    plot(uref(:,2),uref(:,1),'LineWidth',1.5)
    grid on
    xlabel('u velocity')
    ylabel('y')
    legend('SIMPLE','ref','Location','southeast')
    hold off
    axis square;
    
    subplot(1,3,3)
    plot(xc,vb(:,Ny/2),'LineWidth',1.5)
    hold on
    plot(vref(:,1),vref(:,2),'LineWidth',1.5)
    grid on
    xlabel('x')
    ylabel('v velocity')
    legend('SIMPLE','ref')
    hold off
    axis square;

    pause(0.001)
    
end

% plot streamlines :
subplot(1,3,1);
nlines = 5;
startx = linspace(0.62,0.99,nlines);
starty = linspace(0.75,0.975,nlines);

Z0 = get(s,'ZData');
set(s,'ZData',Z0 - 1.1);
stream = streamline(xc,yc,ub',vb',startx,starty,[0.5]);
set(stream,'Color','#EDB120');

% plot 1D cuts:
subplot(1,3,2);
plot(ub(Nx/2,:)',yc,'LineWidth',1.5)
hold on
plot(uref(:,2),uref(:,1),'LineWidth',1.5)
grid on
xlabel('u velocity')
ylabel('y')
legend('SIMPLE','ref','Location','southeast')
axis square;

subplot(1,3,3);
plot(xc,vb(:,Ny/2),'LineWidth',1.5)
hold on
plot(vref(:,1),vref(:,2),'LineWidth',1.5)
grid on
xlabel('x')
ylabel('v velocity')
legend('SIMPLE','ref')
axis square;




