% Semi-implicit scheme of Casulli for the 2D shallow water equations with wetting and drying
clear;
% close all;

global g gamma Nx Ny dx dy dt;

%  physical parameters
g = 9.81;       % acceleration due to gravity
gamma = 0;   % bottom friction coefficient
xL =-5;
xR = 5;
yL =-5;
yR = 5;
% method parameters
Nx = 128;
Ny = 128;
dx = (xR-xL)/Nx;
dy = (yR-yL)/Ny;

x  = linspace(xL,xR,Nx+1); % face x-coordinates of the cells
y  = linspace(yL,yR,Ny+1); % face y-coordinates of the cells
xb = 0.5*(x(2:end) + x(1:end-1)); % cell bary-centers in x 
yb = 0.5*(x(2:end) + y(1:end-1)); % cell bary-centers in y

% time parameters
t    = 0;
tend = 10;
dt   = 0.1;
CFL  = 0.9;    % Courant-Friedrichs-Levi number <= 1

% initial conditions
tmp = 1;
[Xb,Yb] = meshgrid(xb,yb);
[X ,Y ] = meshgrid(x ,y );
etab  = 0.01 + exp(-( Xb.^2 + Yb.^2 )/(2*tmp^2)); % free surface elevation measured from the rest level
bathb = bathymetry(Xb,Yb); % bathymetry in the cell centers
bathx = bathymetry(X ,Yb); % bathymetry on the cell faces in x
bathy = bathymetry(Xb,Y ); % bathymetry on the cell faces in y

u = zeros(Nx+1,Ny); % x-velocity at the cell faces
v = zeros(Nx,Ny+1); % y-velocity at the cell faces 

% plot initial conditions
surf(xb,yb,etab,'EdgeColor','none','FaceColor','interp');
% axis equal;

% t loop
Nt = 100000000;
for n = 1:Nt
    umax = max(max(abs(u)));
    vmax = max(max(abs(v)));
    if (max(umax,vmax) < 10)
        dt = 0.1;
    else
        dt = CFL/( umax/dx + vmax/dy );         % time step restriction of FTCS
    end
    if(t + dt > tend)
        dt = tend - t;
    end
    if(t >= tend)
        break
    end
    if( t >= 3.7)
        zz=1;
    end
    Hb = max(0,etab + bathb);
    % compute total water depth on the cell faces
    [Hmx,Hpx,Hmy,Hpy,rhsb,Hx,Hy] = LinearPartCoef(Hb,u,v); 

    % Newton-type iterations
    kMax = 100;
    for k = 1:kMax
        % determine the wet cells
        Hb = max(0,etab + bathb);
        wet = Hb > 0; %((0.5*(Hpx + Hmx) > 0) + (0.5*(Hpy + Hmy) > 0)) > 0;
        Teta = MatVecProd(etab,bathb,Hmx,Hpx,Hmy,Hpy,wet);
        residual = Teta - rhsb;
        res_norm = norm(residual);
        if(max(max(abs(residual)))<1e-12)
            break
        end
        etab = etab - CGSolver(residual,bathb,Hmx,Hpx,Hmy,Hpy,wet);
    end
    % update the velocity
    [u,v] = VelocityUpdate(u,v,etab,Hx,Hy);

    t = t + dt;
    % plot initial conditions
    surf(xb,yb,etab,'EdgeColor','none','FaceColor','interp');
    axis([xL xR yL yR 0 1])
    title(strcat('t=',num2str(t)))
    pause(0.001)
end

