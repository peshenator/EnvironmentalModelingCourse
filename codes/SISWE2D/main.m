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
Nx = 200;
Ny = Nx;
dx = (xR-xL)/Nx;
dy = (yR-yL)/Ny;

x  = linspace(xL,xR,Nx+1); % face x-coordinates of the cells
y  = linspace(yL,yR,Ny+1); % face y-coordinates of the cells
xb = 0.5*(x(2:end) + x(1:end-1)); % cell bary-centers in x 
yb = 0.5*(x(2:end) + y(1:end-1)); % cell bary-centers in y

% time parameters
t    = 0;
tend = 1;
dt   = 0.1;
CFL  = 0.9;    % Courant-Friedrichs-Levi number <= 1

% initial conditions
tmp = 0.1;
[Xb,Yb] = meshgrid(xb,yb);
[X ,Y ] = meshgrid(x ,y );
% etab  = 0.01 + exp(-( (Xb-0).^2 + (Yb-0).^2 )/(2*tmp^2)); % free surface elevation measured from the rest level

% *** test 5.2 from Dumbser&Casulli, DOI:10.1016/j.amc.2013.02.041
% etab  = 1 + exp(-0.5*( (Xb-0).^2 + (Yb-0).^2 )/tmp^2); % free surface elevation measured from the rest level
% *** test 5.3. Two-dimensional cylindrical dambreak over a bottom step, DOI:10.1016/j.amc.2013.02.041
etab = zeros(Nx,Ny);
Hb   = zeros(Nx,Ny);
bathb= zeros(Nx,Ny);
r0 = 0;
etal = 1.0;
etar = 0.5;
bl = 0.2;
br = 0;
etab = 0.5*(etal + etar) + 0.5*(etar - etal)*erf(Xb'/0.1);
bathb= 0.5*(bl   + br  ) + 0.5*(br   - bl  )*erf(Xb'/0.1);

% etab = 0.5*(etal + etar) + 0.5*(etar - etal)*erf((sqrt(Xb.^2+Yb.^2) - r0)/0.1);
% bathb= 0.5*(bl   + br  ) + 0.5*(br   - bl  )*erf((sqrt(Xb.^2+Yb.^2) - r0)/0.1);
% for i=1:Nx
%     for j=1:Ny
%         if ( sqrt(x(i)^2 + y(j)^2) < r0)
%             Hb(i,j)   = 0.8; 
%             etab(i,j) = etal;
%         else
%             Hb(i,j)   = 0.5;
%             etab(i,j) = etar;
%         end
%     end
% end
% bathb = Hb - etab;
% etab  = etab + 0.01 + exp(-( (Xb+3).^2 + (Yb+4).^2 )/(2*tmp^2)); % free surface elevation measured from the rest level
% bathb = bathymetry(Xb,Yb); % bathymetry in the cell centers

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
    if (max(umax,vmax) < 0.1)
        dt = 0.001;
    else
        dt = CFL/( umax/dx + vmax/dy );         % time step restriction of FTCS
    end
    dt = min(0.01,dt);
    if(t + dt > tend)
        dt = tend - t;
    end
    if(t >= tend)
        break
    end

    Hb = max(0,etab - bathb); % total water depth in the cell centers
    Hub = 0.5*( u(2:Nx+1,:) + u(1:Nx,:) ).*Hb; % x-momentum at the cell centers
    Hvb = 0.5*( v(:,2:Ny+1) + v(:,1:Ny) ).*Hb; % y-momentum at the cell centers
%     [u,v] = MomentumConvectionCons1(Hub,Hvb,Hb,u,v);
%     [u,v] = MCC(Hub,Hvb,Hb,u,v);
    [u,v] = MomentumConvection(u,v);

    % compute total water depth on the cell faces
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
    % plot initial conditions
    subplot(2,1,1)
    surf(xb,yb,(etab)','FaceColor','interp','EdgeColor','none');%'EdgeColor','none',
%     axis([xL xR yL yR 0 2 ])
    clim([min(min(etab)) max(max(etab))])
    xlabel('x')
    ylabel('y')
%     view([0 90])
    subplot(2,1,2)
    plot(xb,etab(:,Ny/2))
    title(strcat('t=',num2str(t)))
    pause(0.001)
end

