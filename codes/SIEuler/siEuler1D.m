3% clear;
% close all;
global Nx dx dt cv rho0 gam uL uR pL pR QL QR plotcall xL xR ;
global MaxPicardP tolP;

% material parameters for Ideal gas Equation of State
cv = 1;
gam = 1.4;
rho0 = 1;

[rhoL, uL, pL, rhoR, uR, pR, xL, xR, tend] = RPinit(3);

% discretization parameters
Nx = 128;
dx = (xR - xL)/Nx;
x  = linspace(xL+dx/2,xR-dx/2,Nx  ); % cell- centers
xv = linspace(xL,xR,Nx+1); % cell edges
CFL = 0.95;

MaxPicardP = 20;
tolP = 1e-5;

% staggered quantities:
rhox = zeros(1,Nx+1);   % rho at the x-faces
ux   = zeros(1,Nx+1);   % u   at the x-faces

%% Initial data and BC: 
rhox(1:Nx/2) = rhoL;    rhox(Nx/2+1:Nx+1) = rhoR;    
ux(1:Nx/2) = uL;        ux(Nx/2+1:Nx+1) = uR;    
mx = rhox.*ux;

% barycenter quantities:
p = zeros(1,Nx);
p(1:Nx/2) = pL;       p(Nx/2+1:Nx) = pR;
rhoE = energy(rhox,ux,p);

QL = [rhoL rhoL*uL rhoE(1) ];
QR = [rhoR rhoR*uR rhoE(Nx)];

time = 0;
plotcall = 0;
auxpar = struct('res',0,'linecolor','#0072BD');
myplot(time,x,xv,rhox,ux,p,rhoE,auxpar);

for n=1:10000000
    v_max =  max(abs(mx./rhox));
%     c_sound = sqrt(gam*pc./(0.5*(rhox(2:Nx+1)+rhox(1:Nx))));
%     v_max = max(v_max + c_sound);
    if (v_max == 0)
        dt = 2e-4;
    else
        dt  = CFL*dx/(v_max);
    end

    if (time >= tend)
        break
    end
    if (time + dt > tend)
        dt = tend  - time;
    end


    %% STEP #1 
    % explicit (convective sub-system)
    %
    % averaging to the cell-centers:
    rho = 0.5*(rhox(2:Nx+1) + rhox(1:Nx));
    m   = 0.5*(  mx(2:Nx+1) +   mx(1:Nx));
    [rho_new,rhox_new] = Convect_q(rho,ux,1);           % convect rho
    [m_new,mx_new]     = Convect_q(  m,ux,2);           % convect momentum m = rho*u
    rhoE_new           = Convect_rhoE(rhoE,rhox,ux,3);  % convect the total energy

    %% STEP #2 
    % implicit (pressure sub-system)
    [rhoE_new,mx_new,residual,iter]=PressureStep(rhox_new,mx_new,rhoE_new);
    
    %% update the solution
    rhox = rhox_new;
    mx   = mx_new;
    ux   = mx_new./rhox_new;

    rhoE = rhoE_new;
    p    = pressure(rhox,ux,rhoE);
    
    time = time + dt;

    %% Plot
    auxpar = struct('res',residual,'linecolor',[0.5 0.0 0.7]);

    subplot(1,3,1)
    % plot(x,rho,'o-','MarkerFaceColor',[1 0 0])
    plot(xv,rhox,'-','Color',auxpar.linecolor)
    grid on
    set(gca,'GridLineStyle',':')
    title({'$\rho$'},'Interpreter','latex')
    
    
    subplot(1,3,2)
    % plot(xv,vv,'o-','MarkerFaceColor',[0 1 0])
    plot(xv,ux,'-','Color',auxpar.linecolor)
    grid on
    set(gca,'GridLineStyle',':')
    title({'$u$'},'Interpreter','latex')
    
    subplot(1,3,3)
    % plot(x,a,'o-','MarkerFaceColor',[0 0.5 0.5])
    plot(x,p,'-','Color',auxpar.linecolor)
    grid on
    set(gca,'GridLineStyle',':')
    title('$p$','Interpreter','latex')
    pause(0.01)

end

Nex = 5*Nx;
[RHO,U,P] = RiemannExact(pL,rhoL,uL,pR,rhoR,uR,tend,xL,xR,Nex,1e-6);
X = linspace(xL,xR,Nex);

subplot(1,3,1)
hold on
% plot(xv,vv,'o-','MarkerFaceColor',[0 1 0])
plot(X,RHO)
hold off

subplot(1,3,2)
hold on
% plot(xv,vv,'o-','MarkerFaceColor',[0 1 0])
plot(X,U)
hold off

subplot(1,3,3)
hold on
% plot(x,a,'o-','MarkerFaceColor',[0 0.5 0.5])
plot(X,P)
hold off






