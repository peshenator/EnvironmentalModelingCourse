% Simplified solidification problem: 
%
% 1) Navier-Stokes equations are discretized using the SIMPLE algorithm
% - viscous terms are implicit  (BTCS scheme)
%
% 2) The nonlinear heat transfer is solved with the BTCS scheme and
% Casulli&Zanolli algorithm of nested Newtonian iterations
%
% 3) The Boussinesq approximation is used for buyoncy forces

clear;
% close all;
% clc;

global Nx Ny dx dy dt nu  uLid g beta TrefBuyoncy lambda uWall vWall TBC xb yb;
global  lambdaS lambdaL hL cS cL rhoL rhoS KS KL Tair Tlake Tc epsilonT
% physical parameters
nu = 1e-2;      % kinematic viscosity
uLid =  0;      % velocity of the lid (top boundary)
uWall = 0;      % u velocity at the walls
vWall = 0;      % v velocity at the walls
g = 9.81;       % accelaration due to gravity
TBC = 2;        % Boundary condition temperature
lambda = 1e-3;  % thermal diffusivity
beta = 1;     % thermal expansion coefficient = drho/dT
Tlake =+0;
Tair  =-1;
TrefBuyoncy = 0;        % reference temperature
epsilonT = 0.05;        % regularization parameter near the solidification temperature


KS = 2.09;     % heat conductivity of the solid phase (ice)
KL = 0.6;      % ... of the liquid phase (water)
hL = 334e-2;         % specific latent heat
rhoS = 0.917;         % density of the solid
rhoL = 1.000;        % density of the liquid
cS   = 210.8;        % heat capacity of the solid
cL   = 418.7;        % heat capacity of the liquid
lambdaS   = KS/(rhoS*cS);
lambdaL   = KL/(rhoL*cL);
Tc   = -0.1;           % critical temperature of the phase change
Pr = nu/lambdaL;     % Prdandtl number for liquid (momentum diff/themaal diff)
Re   = 1/nu;        % Reynolds number
Pe   = Re*Pr;       % Peclet number (advective rate/diffusive rate)

% computational parameters
time = 0;
tend = 10.25;

xL = 0;
xR = 1;
yL = 0;
yR = 1;
Nx = 64;
Ny = Nx;
dx = (xR - xL)/Nx;
dy = (yR - yL)/Ny;

xb = linspace(xL+dx/2,xR-dx/2,Nx);    % x coordinates of the cell-centers
yb = linspace(yL+dy/2,yR-dy/2,Ny);    % y coordinates of the cell-centers
x = linspace(xL,xR,Nx+1);             % x coordinates of the cell-edges
y = linspace(yL,yR,Ny+1);             % y coordinates of the cell-edges
[X,Y] = meshgrid(xb,yb);
CFL = 0.9;                              % CFL number < 1
dt_accel = 1e-3;                        % a number between [0,1], 1 corresponds to the FTCS time step, 0 to fully BTCS

% initial  conditions
u = zeros(Nx+1,Ny);                 % x component of the velocity vector (blue in the lectures) +1 because cell faces
v = zeros(Nx,Ny+1);                 % y component of the velocity vector (green in the lectures)
p = zeros(Nx,Ny);                   % pressure -> cell centre
T = zeros(Nx,Ny);                   % temperature -> cell centre
Tcor= zeros(Nx+1,Ny+1);               % temperature -> cell corners
rhs =zeros(Nx,Ny);                  % right hand-side for the pressure Poisson equation

for i = 1:Nx
    for j = 1:Ny
%        T(i,j) = Tlake - exp(-0.5*( (xb(i) - 0.8 )^2 + (yb(j) - 0.75)^2)/0.1^2 ); %cold bubble
       T(i,j) = TrefBuyoncy + 0.2*exp(-0.5*( (xb(i) - 0.8 )^2 + (yb(j) - 0.2)^2)/0.1^2 ); %hot bubble
    end
end

% TIME LOOP
nmax = 1000000;
for n=1:nmax

    umax = max(max(abs(u)));
    vmax = max(max(abs(v)));

    if max(umax,vmax) < 0.001
        dt = 0.01;
    else
        dt = CFL/( umax/dx + vmax/dy + dt_accel*2*nu*(1/dx^2 + 1/dy^2) );
    end% 
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
    [rhsustar,rhsvstar] = MomConvectionDiffusion(u,v,T);
    % add the grad(pstar) to the velocity
    [rhsustar,rhsvstar] = AddGradP(rhsustar,rhsvstar,rhsustar,rhsvstar,pstar);


    % STEP #3 Implicit:
    Tcor = bary2corners(T,Tlake,Tlake,Tlake,Tair);
    [ustar,erru,ku] = CGsolver(rhsustar,@MatVecProd_u,T,Tcor);
    [vstar,errv,kv] = CGsolver(rhsvstar,@MatVecProd_v,T,Tcor);

    % compute the right hand-side for the Poisson equation
    rhs = (ustar(2:Nx+1,:) - ustar(1:Nx,:))/(dt*dx) + ...
          (vstar(:,2:Ny+1) - vstar(:,1:Ny))/(dt*dy);
    % [pprime,errp,kp] = ConjGradOpt_p(rhs);  
    % use the matrix free Conjugate Gradient method to solve the pressure Poisson equation
    [pprime,errp,kp] = CGsolver(rhs,@MatVecProd_p,T,Tcor);

    % STEP #4 Div-free velocity correction (add grad(Pprime) to the velocity)
    [u,v] = AddGradP(u,v,ustar,vstar,pprime);
    
    % update pressure adding the correction P^{n+1} = P* + P' :
    p = pstar + pprime;
    
    % Temperature convection-diffusion
    T = TConvection(u,v,T);
    % Tempreature nonlinear diffusion
    T = TCasulliZanolli(T);
    
    % time update
    time = time + dt;
    
    % velocity at the cell barybenters (only for the vizualization!)
    divuv = (u(2:Nx+1,:) - u(1:Nx,:))/dx + (v(:,2:Ny+1) - v(:,1:Ny))/dy;
    ub = 0.5*( u(2:Nx+1,:) + u(1:Nx,:) );
    vb = 0.5*( v(:,2:Ny+1) + v(:,1:Ny) );
    
    subplot(1,2,1)
    hold off
    surf(xb,yb,sqrt(ub.^2+vb.^2)','EdgeColor','none','FaceColor','interp'); %sqrt(ub.^2+vb.^2)
    view([0 90]) % position the camera, vision angle
    hold on 
    title(strcat('Vertical velocity, time = ',num2str(time)))
    xlabel('x [m]')
    ylabel('y [m]')
    axis square;
    colorbar;

    subplot(1,2,2)
    s = surf(xb,yb,T','EdgeColor','none','FaceColor','interp'); %sqrt(ub.^2+vb.^2)
    view([0 90]) % position the camera, vision angle
    hold on
    contour(X,Y,T',[ Tc Tc ],'b', 'LineWidth', 2)
    quiver(xb(1:2:Nx),yb(1:2:Ny),ub(1:2:Nx,1:2:Ny)',vb(1:2:Nx,1:2:Ny)','k');
    hold off

    title(strcat('Temperature and velocity field, time = ',num2str(time)))
    xlabel('x [m]')
    ylabel('y [m]')
    % clim([-50 0])
    Z0 = get(s,'ZData');
    set(s,'ZData',Z0 - 0.3)
    axis square;
    colorbar;

    pause(0.01)
    
end




