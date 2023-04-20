% Semi-implicit scheme of Casulli for the 2D shallow water equations with wetting and drying
clear;
% close all;

global g gamma Nx Ny dx dy dt;

%  physical parameters
xL =-5; xR = 5; 
yL =-5; yR = 5;
g = 9.81;       % acceleration due to gravity
gamma = 0;   % bottom friction coefficient

% method parameters
Nx = 128;
Ny = Nx;
dx = (xR-xL)/Nx;
dy = (yR-yL)/Ny;

x  = linspace(xL,xR,Nx+1); % face x-coordinates of the cells
y  = linspace(yL,yR,Ny+1); % face y-coordinates of the cells
xb = 0.5*(x(2:end) + x(1:end-1)); % cell bary-centers in x 
yb = 0.5*(x(2:end) + y(1:end-1)); % cell bary-centers in y

% time parameters
t   = 0;
dt  = 0.1;
CFL = 0.9;    % Courant-Friedrichs-Levi number <= 1

% initial conditions
[Xb,Yb] = meshgrid(xb,yb);
[X ,Y ] = meshgrid(x ,y );
% *** test #1: Smoothe wave (test 5.2 from Dumbser&Casulli, DOI:10.1016/j.amc.2013.02.041)
% tend = 0.15;
% dt0  = 0.002;
% etab  = 1 + exp(-0.5*( (Xb-0).^2 + (Yb-0).^2 )/0.1^2); % free surface elevation measured from the rest level
% bathb = zeros(Nx,Ny);

% *** test #2: Smooth wave with wetting/drying
% dt0 = 0.01;
% tend = 10;
% etab  = 0.01 + exp(-( (Xb-0).^2 + (Yb-0).^2 )/(2*0.1^2)); % free surface elevation measured from the rest level
% bathb = zeros(Nx,Ny);

% *** test #3: 1D dam break over a bottom with a step
% tend = 1;
% dt0  = 0.01;
% etal = 1.0;
% etar = 0.5;
% bl   = 0.2;
% br   = 0;
% etab = 0.5*(etal + etar) + 0.5*(etar - etal)*erf(Xb'/0.1);
% bathb= 0.5*(bl   + br  ) + 0.5*(br   - bl  )*erf(Xb'/0.1);

% *** test #4: Two-dimensional cylindrical dambreak over a bottom step
% tend = 0.2;
% dt0=2e-3;
% r0   = 1;
% etal = 1.0;
% etar = 0.5;
% bl   = 0.2;
% br   = 0;
% etab = 0.5*(etal + etar) + 0.5*(etar - etal)*erf((sqrt(Xb.^2+Yb.^2) - r0)/0.1);
% bathb= 0.5*(bl   + br  ) + 0.5*(br   - bl  )*erf((sqrt(Xb.^2+Yb.^2) - r0)/0.1);

% *** test #5:
dt0 = 0.01;
tend = 10;
etab = 0.01 + 1.5*exp(-( (Xb-0).^2 + (Yb-0).^2 )/(2*0.4^2));%0.1*Xb';
bathb = -1 + 0.1*Xb'.^2 + 0.1*Yb'.^2;

u = zeros(Nx+1,Ny); % x-velocity at the cell faces
v = zeros(Nx,Ny+1); % y-velocity at the cell faces 

% plot initial conditions
surf(xb,yb,etab','EdgeColor','none','FaceColor','interp');
% axis equal;

% t loop
Nt = 100000000;
for n = 1:Nt
    umax = max(max(abs(u)));
    vmax = max(max(abs(v)));
    dt = CFL/( umax/dx + vmax/dy );         % time step restriction of FTCS
    dt = min(dt0,dt);
    if(t + dt > tend)
        dt = tend - t;
    end
    if(t >= tend)
        break
    end

    % Explicit subsystem (velocity convection)
    [u,v] = MomentumConvection(u,v);

    Hb = max(0,etab - bathb); % total water depth in the cell centers
    [Hmx,Hpx,Hmy,Hpy,rhsb,Hx,Hy] = LinearPartCoef(Hb,u,v); 

    % Newton-type iterations
    kMax = 100;
    for k = 1:kMax
        % determine the wet cells
        Hb = max(0,etab - bathb);
        wet = Hb > 0;
        Teta = Hb + MatVecProd(etab,bathb,Hmx,Hpx,Hmy,Hpy,wet);
        residual = Teta - rhsb;
        res_norm = norm(residual);
        if(max(max(abs(residual)))<1e-10)
            break
        end
        [deta,CGk,CGerr] = CGSolver(residual,bathb,Hmx,Hpx,Hmy,Hpy,wet);
        etab = etab - deta;
    end
    % update the velocity
    [u,v] = VelocityUpdate(u,v,etab,Hx,Hy);

    t = t + dt;
    
    %% Postprocessing and plotting
    etabPlot = etab;
    etabPlot(Hb<1e-2)=NaN;

    subplot(2,1,1)
    sb = surf(xb,yb,bathb','EdgeColor','none','FaceColor','flat');
    alpha 0.25;
    hold on
    surf(xb,yb,(etabPlot)','FaceColor','interp','EdgeColor','none');%'EdgeColor','none',
    camlight;
    axis([xL xR yL yR -1 0.5])
    clim([min(min(etab)) max(max(etab))])
    xlabel('x')
    ylabel('y')
    hold off

%     view([0 90])
    subplot(2,1,2)
    plot(xb,etabPlot(:,Ny/2))
    hold on
    plot(xb,bathb(:,Ny/2),'r')
    title(strcat('t=',num2str(t)))
    hold off
    pause(0.001)
end

