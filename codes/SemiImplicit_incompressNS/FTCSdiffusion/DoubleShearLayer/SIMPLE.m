% SIMPLE scheme for the incompressible Navier-Stokes equations
clear;
% close all;
global imax jmax dt dx dy nu g beta uLid T0 TB lambda xb yb;

set(0, 'DefaultLineLineWidth', 1.25,'defaultAxesXGrid','on','defaultAxesYGrid','on');

% model parameters
nu     = 1e-5;     % kinematic viscosity
uLid   = 1;        % x-velocity of the lid
g      = 9.8;      % acceleration due to the gravity
beta   = 0;        % beta = dT/drho in the Boussinesq approximation
T0     = 0;        % reference temperature
TB     = 0;        % temperature at the boundary
lambda = 1e-3;  % heat conduction coefficient

% computational domain
xL = 0;
xR = 1;
yL = 0;
yR = 1;

imax = 128;
jmax = 128;

dx = (xR - xL)/imax;
dy = (yR - yL)/jmax;

x = linspace(xL,xR,imax+1);     % coordinates of the cell edges
y = linspace(yL,yR,jmax+1);     % coordinates of the cell edges

% xb = linspace(xL,xR,imax);    % coordinates of the barycenters of the cells
% yb = linspace(yL,yR,jmax);    % coordinates of the barycenters of the cells

xb = linspace(xL+dx/2,xR-dx/2,imax);    % coordinates of the barycenters of the cells
yb = linspace(yL+dy/2,yR-dy/2,jmax);    % coordinates of the barycenters of the cells


t = 0;
tend = 3;
CFL = 0.75;      % CFL number for the explicit step (step number 2)

% initial coditions
u = zeros(imax+1,jmax);     % x-component of the velocity (blue in the lecture)
v = zeros(imax,jmax+1);     % y-component of the velocity (green in the lecture)

ub = zeros(imax,jmax);      % barycenter u-velocty (only for ploting)
vb = zeros(imax,jmax);      % barycenter v-velocty (only for ploting)

P = zeros(imax,jmax);       % the pressure in the barycenters (red)

% 
rhs = zeros(imax,jmax);     % right-hand side for the pressure Poisson eq.
T   = zeros(imax,jmax);     % Temperature in the barycenters
Id = eye(imax*jmax);
%% Initial data from Michael
u1 = 0.5;
u2 =-0.5;
um = (u1-u2)/2;
L = 0.025;
% u
for j = 1:jmax
    if (yb(j) <= 0.25)
        u(:,j) = u1 - um*exp(( yb(j)-0.25 )/L);
    elseif(yb(j)<0.5 && yb(j) >= 0.25)
        u(:,j) = u2 + um*exp((-yb(j)+0.25 )/L);
    elseif(yb(j)<0.75 && yb(j) >= 0.5)
        u(:,j) = u2 + um*exp(( yb(j)-0.75 )/L);
    elseif(yb(j)<1.0 && yb(j) >= 0.75)
        u(:,j) = u1 - um*exp((-yb(j)+0.75 )/L);
    end
end
% v
for j = 1:jmax+1
    v(:,j) = 1e-2*sin(4*pi*xb(:));
end
uc = fx2c(u);
vc = fy2c(v);
u = c2fx(uc);
v = c2fy(vc);

%% Initial data from Walter
% theta = 30;
% delta = 0.05;
% for j = 1:jmax
%     if (yb(j) <= 0.5)
%         u(:,j) = tanh(theta*(yb(j) - 0.25));
%     else
%         u(:,j) = tanh(theta*(0.75 - yb(j)));
%     end
% end
% for j = 1:jmax+1
%     v(:,j) = delta*sin(2*pi*xb(:));
% end
% 
% for i = 1:imax
%     if (xb(i) <= 0.5)
%         v(i,:) = 0.5*tanh(theta*(xb(i) - 0.25));
%     else
%         v(i,:) = 0.5*tanh(theta*(0.75 - xb(i)));
%     end
% end
% for i = 1:imax+1
%     u(i,:) = delta*sin(2*pi*yb(:));
% end

%%

% fig1 = figure;
% set(gcf,'Position',[100 100 1050 500])
    uc = 0.5*( u(2:imax+1,:) + u(1:imax,:) );
    vc = 0.5*( v(:,2:jmax+1) + v(:,1:jmax) );
    subplot(1,2,1)
    s = surf(xb,yb,uc','EdgeColor','none','FaceColor','interp');
%     view([0 90]) % position the camera, vision angle
    hold on 
%     quiver(xb,yb,uc',vc','k');
    title(strcat('Current time = ',num2str(t)))
    xlabel('x [m]')
    ylabel('y [m]')
%     Z0 = get(s,'ZData');
%     set(s,'ZData',Z0 - 10)
    axis square;
    hold off

    subplot(1,2,2)
    s = surf(xb,yb,vc','EdgeColor','none','FaceColor','interp');
%     view([0 90]) % position the camera, vision angle
    hold on 
%     quiver(xb,yb,uc',vc','k');
    title(strcat('Current time = ',num2str(t)))
    xlabel('x [m]')
    ylabel('y [m]')
%     Z0 = get(s,'ZData');
%     set(s,'ZData',Z0 - 10)
    axis square;
    hold off
     
% time loop
for n=1:100000000
    % compute the time step dt:
    umax = max( max( abs(u) ) );
    vmax = max( max( abs(v) ) );
    if max(umax,vmax) < 1e-5
        dt = 1e-3;
    else
        dt = CFL/( umax/dx + vmax/dy + 2*nu*(1/dx^2 + 1/dy^2) );
    end
    
    if (t + dt > tend)
        dt = tend - t;
    end
    if (t >= tend)
        break
    end

    % *** STEP 1 ***
    Pstar = P;      % initial guess for the pressure

    % *** STEP 2 ***
    [ustar,vstar] = MomentumConvectionDiffusion(u,v,T);         % update ustar, vstar (convection + diffusion)
    [ustar,vstar] = MomentumPressureGrad(ustar,vstar,Pstar);    % update ustar,vstar (add the pressure gradient)
%     Tnew = TemperatureConvectionDiffusion(u,v,T);               % update Temperature (convection + diffusion)

    % *** STEP 3 ***
    rhs = DivVelocityStar(ustar,vstar);     % compute divergence of the (ustar,vstar)
    opts.SYM=true;
%     opts.POSDEF=true;
%     [PM, RHS] = PoissonMatrixPeriodicBC(Pstar,rhs,imax,jmax,dx^2,dy^2);
%     X = linsolve(1e-7*Id + PM,RHS,opts);
%     Pprime1 = reshape(X,imax,jmax);


    [Pprime,CGerr,CGiter] = CG_matrix_free(rhs);           % solve the Poisson equation for the Pprime

    % *** STEP 4 ***
%     [u,v] = MomentumPressureGrad(ustar,vstar,Pprime);
    [u,v] = DivFreeVelocityCorrection(u,v,ustar,vstar,Pprime);
    divuv = (u(2:imax+1,:) - u(1:imax,:))/dx + (v(:,2:jmax+1) - v(:,1:jmax))/dy;

    P = Pstar + Pprime;     % update the pressure
%     T = Tnew;               % update the temperature
   
    t = t + dt;   % advance time
%% Plot part
    % vorticity (only z-component is non-zero) 
    omegaz = vorticity(u,v,imax,jmax,dx,dy);
    % for plotting only, average the velocities into the barybenter points
    ub = 0.5*(u(2:imax+1,:)+u(1:imax,:));
    vb = 0.5*(v(:,2:jmax+1)+v(:,1:jmax));
    hold off
    %
    subplot(1,2,1)
%     hold off
    s1 = surf(x,y,omegaz','EdgeColor','none','FaceColor','interp');
    view([0 90]) % position the camera, vision angle
%     hold on 
%     quiver(xb,yb,uc',vc','k');
    title(strcat('Current time = ',num2str(t)))
    xlabel('x [m]')
    ylabel('y [m]')
%     Z0 = get(s,'ZData');
%     set(s,'ZData',Z0 - 10)
    axis square;

    subplot(1,2,2)
%     hold off
    s2 = surf(xb,yb,divuv','EdgeColor','none','FaceColor','interp');
    view([0 90]) % position the camera, vision angle
%     hold on 
%     quiver(xb,yb,uc',vc','k');
    title(strcat('Current time = ',num2str(t)))
    xlabel('x [m]')
    ylabel('y [m]')
%     Z0 = get(s,'ZData');
%     set(s,'ZData',Z0 - 10)
    axis square;

   
    pause(0.001)
end

% plot streamlines :
nlines = 5;
startx = linspace(0.62,0.99,nlines);
starty = linspace(0.75,0.975,nlines);

% fig2 = figure('Position',[16,268,399,300]);
% view(2)
% s.EdgeColor = 'none';
% colorbar;
% hold on
Z0 = get(s,'ZData');
set(s,'ZData',Z0 - 1.1);
stream = streamline(xb,yb,ub',vb',startx,starty,[0.5]);
set(stream,'Color','#EDB120');

% plot 1D cuts:
fig2 = figure('Position',[926 483 909 369]);
subplot(1,2,1)
plot(ub(imax/2,:)',yb)
grid on
xlabel('u velocity')
ylabel('y')

subplot(1,2,2)
plot(xb,vb(:,jmax/2))
grid on
xlabel('x')
ylabel('v velocity')



