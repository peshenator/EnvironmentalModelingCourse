% SIMPLE scheme for the incompressible Navier-Stokes equations
clear;
% close all;
global Nx Ny dt dx dy nu g beta uLid T0 TB lambda xb yb;

set(0, 'DefaultLineLineWidth', 1.25,'defaultAxesXGrid','on','defaultAxesYGrid','on');

% model parameters
nu     = 2e-4;     % kinematic viscosity
uLid   = 0;        % x-velocity of the lid
g      = 9.8;      % acceleration due to the gravity
beta   = 0;        % beta = dT/drho in the Boussinesq approximation
T0     = 0;        % reference temperature
TB     = 0;        % temperature at the boundary
lambda = nu;  % heat conduction coefficient

% computational domain
xL = 0;
xR = 1;
yL = 0;
yR = 1;

Nx = 128;
Ny = Nx;

dx = (xR - xL)/Nx;
dy = (yR - yL)/Ny;

x = linspace(xL,xR,Nx+1);     % coordinates of the cell edges
y = linspace(yL,yR,Ny+1);     % coordinates of the cell edges

xb = linspace(xL+dx/2,xR-dx/2,Nx);    % coordinates of the barycenters of the cells
yb = linspace(yL+dy/2,yR-dy/2,Ny);    % coordinates of the barycenters of the cells


t = 0;
tend = 1.8;
CFL = 0.95;      % CFL number for the explicit step (step number 2)

% initial coditions
u = zeros(Nx+1,Ny);     % x-component of the velocity (blue in the lecture)
v = zeros(Nx,Ny+1);     % y-component of the velocity (green in the lecture)

ub = zeros(Nx,Ny);      % barycenter u-velocty (only for ploting)
vb = zeros(Nx,Ny);      % barycenter v-velocty (only for ploting)

P = zeros(Nx,Ny);       % the pressure in the barycenters (red)
T = zeros(Nx,Ny);

% 
rhs = zeros(Nx,Ny);     % right-hand side for the pressure Poisson eq.

%% Initial data for the double shear layer
[u,v,T] = DoubleShearLayer2(u,v,T,1);
% [u,v,T] = BubbleAdvection(u,v,T,[1 1]);
% Tex = T;

% fig1 = figure;
% set(gcf,'Position',[100 100 1050 500])
    uc = 0.5*( u(2:Nx+1,:) + u(1:Nx,:) );
    vc = 0.5*( v(:,2:Ny+1) + v(:,1:Ny) );
    subplot(1,2,1)
    s = surf(xb,yb,uc','EdgeColor','none','FaceColor','interp');
    view([0 90]) % position the camera, vision angle
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
    s = surf(xb,yb,T','EdgeColor','none','FaceColor','interp');
    view([0 90]) % position the camera, vision angle
    colormap(gca,jet)
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
    [ustar,vstar] = MomentumConvectionDiffusion(u,v,T);      % update ustar, vstar (convection + diffusion)
    [ustar,vstar] = MomentumPressureGrad(ustar,vstar,Pstar); % update ustar,vstar (add the pressure gradient)
    
    % *** STEP 3 ***
    rhs = DivVelocityStar(ustar,vstar);     % compute divergence of the (ustar,vstar)
    [Pprime,CGerr,CGiter] = CG_matrix_free(rhs);           % solve the Poisson equation for the Pprime
%     Pprime = Pstar;

    % *** STEP 4 ***
    [u,v] = DivFreeVelocityCorrection(u,v,ustar,vstar,Pprime);
    divuv = (u(2:Nx+1,:) - u(1:Nx,:))/dx + (v(:,2:Ny+1) - v(:,1:Ny))/dy;

    P = Pstar + Pprime;     % update the pressure
    T = ScalConvectionDiffusionMUSCL(u,v,T,dt);               % update Temperature (convection + diffusion)

   
    t = t + dt;   % advance time
%% Plot part
    % vorticity (only z-component is non-zero) 
    omegaz = vorticity(u,v,Nx,Ny,dx,dy);
    % for plotting only, average the velocities into the barybenter points
    ub = 0.5*(u(2:Nx+1,:)+u(1:Nx,:));
    vb = 0.5*(v(:,2:Ny+1)+v(:,1:Ny));
    hold off
    %
    subplot(1,2,1)
%     hold off
    surf(x,y,omegaz','EdgeColor','none','FaceColor','interp');
%     contour(x,y,-omegaz',30);
    view([0 90]) % position the camera, vision angle
    colormap(gca,"pink")
    caxis([-25 25])
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
    surf(xb,yb,T','EdgeColor','none','FaceColor','interp');
    view([0 90]) % position the camera, vision angle
    colormap(gca,"copper")
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



