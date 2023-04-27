% Semi-implicit scheme of Casulli for the 2D shallow water equations with
% wetting and drying and sediment transport
clear;
% close all;

global g gamma Nx Ny dx dy dt;
global sA sm sUc sphi;

%  physical parameters
xL =-5; xR = 5; 
yL =-5; yR = 5;
g = 9.81;       % acceleration due to gravity
gamma = 0.0;   % bottom friction coefficient

sA   = 0.24;  % some coefficient in the sediment transport flux
sm   = 3;    % sediment flux parameter
sUc  = 0.3;  % critical velocity below which sediment does not move
sphi = 0.05; % porosity of the sediment

% method parameters
Nx = 64;
Ny = Nx;
dx = (xR-xL)/Nx;
dy = (yR-yL)/Ny;

x  = linspace(xL,xR,Nx+1); % face x-coordinates of the cells
y  = linspace(yL,yR,Ny+1); % face y-coordinates of the cells
xb = 0.5*(x(2:end) + x(1:end-1)); % cell bary-centers in x 
yb = 0.5*(x(2:end) + y(1:end-1)); % cell bary-centers in y

% time parameters
t   = 0;
dt0 = 0.1;
CFL = 0.5;    % Courant-Friedrichs-Levi number <= 1

% initial conditions
[Xb,Yb] = meshgrid(xb,yb);
[X ,Y ] = meshgrid(x ,y );
dt0  = 0.1;
tend = 100;
etab = 1*ones(Nx,Ny) + 0*exp(-0.5*( (Xb-0).^2 + (Yb-0).^2 )/0.9^2);
bathb= 0.4*exp(-0.5*( (Xb-0).^2 + (Yb-0).^2 )/0.9^2);
b0 = bathb;
Hb = max(0,etab - bathb); % total water depth in the cell centers
u = 1.0*ones(Nx+1,Ny); % y-velocity at the cell faces 
v = zeros(Nx,Ny+1); % y-velocity at the cell faces 

% plot initial conditions
subplot(2,1,1)
surf(xb,yb,(etab)','FaceColor','interp','EdgeColor','none');%'EdgeColor','none',
alpha 0.25;
hold on
sb = surf(xb,yb,bathb','EdgeColor','none','FaceColor','flat');
camlight;
axis([xL xR yL yR -1 1.5])
xlabel('x')
ylabel('y')
hold off
% axis equal;

% t loop
Nt = 100000000;
for n = 1:Nt
    umax = max(max(abs(u)));
    vmax = max(max(abs(v)));
    bmax = max(max(abs(Diffbflux(u,v,Hb))));
    Cmax = max([ umax vmax bmax]);
    dt = min(dt0,CFL*dx/Cmax)
    if(t + dt > tend)
        dt = tend - t;
    end
    if(t >= tend)
        break
    end

    %
    Hb = max(0,etab - bathb); % total water depth in the cell centers
    [Hmx,Hpx,Hmy,Hpy,Hx,Hy,tHx,tHy,rhsb] = LinearPartCoef(Hb,u,v); 

    % Explicit subsystem (velocity convection)
    [u,v] = MomentumConvection(u,v);
    [u,v] = VelocityFilter(u,v,tHx,tHy);
    % bathymetry update
    bathb = BathUpdate(bathb,u,v,Hb);
  
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
    epsPlot = 1e-2;
    etabPlot = etab;
    etabPlot(Hb<epsPlot)=NaN;
    ub = 0.5*(u(1:end-1,:) + u(2:end,:));
    vb = 0.5*(v(:,1:end-1) + v(:,2:end));
    uvb = sqrt(ub.^2 + vb.^2);
    uvb(Hb<epsPlot)=NaN;

    subplot(2,2,1)
    surf(xb,yb,(etabPlot)','FaceColor','interp','EdgeColor','none');%'EdgeColor','none',
    alpha 0.75;
    hold on
    surf(xb,yb,bathb','EdgeColor','none','FaceColor','flat');
    camlight;
    axis([xL xR yL yR -0.05 1.1])
    clim([min(min(etab))-1e-3 max(max(etab))+1e-3])
    xlabel('x')
    ylabel('y')
    zlabel('$b,\ \eta$','Interpreter','latex')
    hold off
    frame(n) = getframe;

%     view([0 90])
    subplot(2,2,2)
    plot(xb,etabPlot(:,Ny/2),'o')
    hold on
    plot(xb,ub(:,Ny/2),'*')
    plot(xb,b0(:,Ny/2))
    plot(xb,bathb(:,Ny/2),'r')
    xlabel('x')
    ylabel('$b,\ \eta,\ u$','Interpreter','latex')
    title(strcat('t = ',num2str(t),',   dt = ',num2str(dt)))
    hold off

    subplot(2,2,3)
    plot(xb,etabPlot(:,Ny/2),'o')
    axis([xL xR .50 1.05])
    xlabel('x')
    ylabel('\eta')
    pause(0.001)
end