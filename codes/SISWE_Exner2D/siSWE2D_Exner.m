% Semi-Implicit method of Casulli (see the paper in the PDF file) for the
% SWE with wetting and drying
 clear;

global Nx Ny dt dx dy gamma g uWall;
global sA sm sUc sphi;

% physical parameters
 xL =-5;    xR = 5;
 yL =-2.5;  yR = 2.5;
 g = 9.81;      % acceleration due to gravity
 gamma = 0;     % bottom friction = const in this code

sA   = 0.14;  % some coefficient in the sediment transport flux
sm   = 3;    % sediment flux parameter
sUc  = 0.1;  % critical velocity below which sediment does not move
sphi = 0.05; % porosity of the sediment

% numerical parameters
Nx = 80; Ny = 40;
dx = (xR - xL)/Nx;  dy = (yR - yL)/Ny;
x = linspace(xL,xR,Nx+1);
y = linspace(yL,yR,Ny+1);
xb = linspace(xL+dx/2,xR-dx/2,Nx);
yb = linspace(yL+dy/2,yR-dy/2,Ny);

% time parameters
t   = 0;
dt  = 0.1;
CFL = 0.5;    % Courant-Friedrichs-Levi number <= 1

% initial conditions
[Xb,Yb] = meshgrid(xb,yb);
[X ,Y ] = meshgrid(x ,y );

tend = 100;
dt0  = 0.1;
etab  = ones(Nx,Ny); % free surface elevation measured from the rest level
% set initial bathymetry
bathb=  0.4*exp(-0.5*( (Xb'+0).^2 + (Yb'-0).^2 )/0.5^2)...
      + 0.8*exp(-0.5*( (Xb'-2.5).^2 + (Yb'+1).^2 )/0.5^2);
bathb0 = bathb;      % remember initial bathymetry for comparison

uWall = 0.5;
u = uWall*ones(Nx+1,Ny);
v = zeros(Nx,Ny+1);
ub = 0.5*( u(2:Nx+1,:) + u(1:Nx,:) );   % u at the barycenters
vb = 0.5*( v(:,2:Ny+1) + v(:,1:Ny) );   % v at the barycenters

Hb = max(0,etab-bathb);
Hb0 = sum(sum(Hb))/(Nx*Ny);
Hu0 = 0;

% time loop
for n = 1:1000000
    umax = max(max(abs(u)));
    vmax = max(max(abs(v)));
    cxmax = max(max(max(abs(ub - sqrt(g*Hb)),abs(ub + sqrt(g*Hb)))));    % max of the characteristic speed for the full flux
    cymax = max(max(max(abs(vb - sqrt(g*Hb)),abs(vb + sqrt(g*Hb)))));    % max of the characteristic speed for the full flux
    [aBx,aBy] = aBath(u,v,Hb);
    aBxmax = max(max(abs(aBx)));
    aBymax = max(max(abs(aBy)));
    Fr  = max(umax,vmax)/max(max(Hb));     % Froud number (equivalent of the Mach number for the Navier-Stokes(Euler) equations) 
    dte = CFL/( cxmax/dx + cymax/dy );     % time step for a fully EXPLICIT scheme
    dt  = CFL/(  umax/dx + vmax/dy  );     % time step for the semi-implicit scheme
    dt  = min(dt0,dt);
    if(t + dt > tend)
        dt = tend - t;
    end
    if(t >= tend)
        break
    end

    % compute the deapth H
    Hb = max(0,etab - bathb);
    [Hxm,Hxp,Hym,Hyp,Hx,Hy,tHx,tHy,rhs] = LinearPartCoeff(Hb,u,v);

    % *** Step 1 (Explicit)
    [u,v] = MomentumConvection(u,v);    % compute u* and v*, see the lecture notes notations 
    [u,v] = VelocityFilter(u,v,Hb);     % apply a filter to velocities in the dry cells
    
    bathb = BathymetryUpdate(bathb,u,v,Hb);

    % *** Step 2 (Implicit), update eta and then the velocity
    kMax = 100;
    for k = 1:kMax
        Hb = max(0,etab-bathb);
        wet = Hb > 0;   % this is what is called porosity p(x,y,z) in the lecture
        Meta = Hb + MatVecProd(etab,bathb,Hxm,Hxp,Hym,Hyp,wet);
        residual = Meta - rhs;  % rhs = b in the lecture
        residual_norm = max(max(abs(residual)));
        if ( residual_norm < 1e-8 )
            break;
        end
        [detab,CGerr,CGk] = CGsolver(residual,@MatVecProdNewton,bathb,Hxm,Hxp,Hym,Hyp,wet);
        etab = etab - detab;
    end
    % Now, updtae the velocity
    [u,v] = VelocityUpdate(u,v,etab,Hx,Hy);
    % update time
    t = t + dt;
    
    % postprocessing part
    epsPlot  = 1e-10;
    etabPlot = etab;
    etabPlot( Hb < epsPlot ) = NaN;
    ub = 0.5*( u(2:Nx+1,:) + u(1:Nx,:) );
    vb = 0.5*( v(:,2:Ny+1) + v(:,1:Ny) );
    Herror = abs(Hb0 - sum(sum(Hb))/(Nx*Ny) )/Hb0;
    Hu_error = abs(Hu0 - sum(sum( Hb.*ub ))/(Nx*Ny));


    % plotting part
    subplot(121)
    s = surf(xb,yb,etabPlot','EdgeColor','none','FaceColor','flat');
    colormap("hsv")
    alpha 0.5;
    view([-25 10])
    hold on
    surf(xb,yb,bathb','FaceColor','interp','EdgeColor','none')
    camlight;
    axis( [xL xR yL yR -1  1.5 ] )
    clim( [ min(min(etab)) - 1e-3,max(max(etab)) + 1e-3 ]  );
    xlabel('x')
    ylabel('y')
    hold off
    title(strcat('H_{error} = ',num2str(Herror)))
    frame(n) = getframe;

    subplot(122)
    [y1,j1] = min(abs(yb + 1));
    [y2,j2] = min(abs(yb - 1));
    plot(xb,etabPlot(:,j1),'o')
    hold on
    plot(xb,etabPlot(:,j2),'o')
    plot(xb,bathb(:,j1))
    plot(xb,bathb(:,j2))
    plot(xb,bathb0(:,j1),'--')
    plot(xb,bathb0(:,j2),'--')
    title(strcat('$t = $',num2str(t),', \ Fr = ',num2str(Fr),'$, \ \Delta t_{imp}/\Delta t_{exp} =$',num2str(dt/dte)),'Interpreter','latex')
    hold off
    legend('\eta(y=-1)','\eta(y=+1)','b(y=-1)','b(y=+1)','b(y=-1,t=0)','b(y=+1,t=0)','Location','west')

    pause(0.001)

end


% Tasks
% 1) check momentum consrvation
% 2) Play with the filter for H = 0







