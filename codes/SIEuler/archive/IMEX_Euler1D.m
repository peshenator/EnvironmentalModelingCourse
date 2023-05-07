clear;
% close all;
global Nx dx dt cv rho0 gam uL uR pL pR QL QR plotcall xL xR ;
global MaxPicardP tolPicard numFlux;

% material parameters for Ideal gas Equation of State
cv = 1;
gam = 1.4;
rho0 = 1;

[rhoL, uL, pL, rhoR, uR, pR, xL, xR, tend] = RPinit(1);

% discretization parameters
numFlux = 1;    % 1 = Rusanov, 2 = FORCE
Nx = 1024;
dx = (xR - xL)/Nx;
xb = linspace(xL+dx/2,xR-dx/2,Nx  ); % cell- centers
x  = linspace(xL,xR,Nx+1); % cell edges
CFL = 0.1;
beta = 1 - sqrt(2)/2;


MaxPicardP = 20;
tolPicard  = 1e-8;

%% Initial data and BC: 

% cell x-face quantities:
qx = zeros(3,Nx+1);

qx(1,1:Nx/2) = rhoL;        qx(1,Nx/2+1:Nx+1) = rhoR;    
qx(2,1:Nx/2) = rhoL*uL;     qx(2,Nx/2+1:Nx+1) = rhoR*uR;
px(  1:Nx/2) = pL;          px(  Nx/2+1:Nx+1) = pR;
qx(3,:) = rhoE_x(qx,px);

% barycenter quantities:
qb = 0.5*(qx(:,2:Nx+1) + qx(:,1:Nx));
pb = 0.5*(px(2:Nx+1)   + px(1:Nx));

QL = [rhoL rhoL*uL qx(3,1   )]';
QR = [rhoR rhoR*uR qx(3,Nx+1)]';

time = 0;
plotcall = 0;
auxpar = struct('res',0,'linecolor','#4DBEEE');
myplot(time,xb,x,qx(1,:),qx(2,:)./qx(1,:),pb,qb(3,:),auxpar);

for n=1:10000000
    v_max =  max(abs(qx(2,:)./qx(1,:)));
    c_sound = sqrt(gam*pb./qb(1,:));
    c_max = max(v_max + c_sound);
    dtc  = CFL*dx/(c_max);
    if (n < 50)
        dt = 2e-5;
    else
        dt  = CFL*dx/(v_max);
    end
%     dt = min(1e-3,dt);

    if (time >= tend)
        break
    end
    if (time + dt > tend)
        dt = tend  - time;
    end
%     beta = 1;
    [Fqb,Fqx] = ExplUpdate(qb,qx,qb,qx,beta*dt,dx);   % convect all quantities at the cell-centeres
    [qb_star,qx_star,residual,iter] = PressureSubSysCG(Fqb,Fqx,beta*dt,dx);
%     qb = qb_star;
%     qx = qx_star;

    [Fqb,Fqx] = ExplUpdate(qb,qx,qb,qx,(beta-1)*dt,dx);
    [Fqb,Fqx] = ExplUpdate(Fqb,Fqx,qb_star,qx_star,(2-beta)*dt,dx);
    [Fqb,Fqx] = ImplUpdate(Fqb,Fqx,qb_star,qx_star,(1-beta)*dt,dx);
    [qb,qx,residual,iter] = PressureSubSysCG(Fqb,Fqx,beta*dt,dx);
    time = time + dt
    
    %% Plot
    pb = pressure(qb(1,:),qb(2,:)./qb(1,:),qb(3,:));

    subplot(1,3,1)
    plot(x,qx(1,:),'-','Color',auxpar.linecolor)
    grid on
    set(gca,'GridLineStyle',':')
    title({'$\rho$'},'Interpreter','latex')
    
    
    subplot(1,3,2)
    plot(x,qx(2,:)./qx(1,:),'-','Color',auxpar.linecolor)
    grid on
    set(gca,'GridLineStyle',':')
    title({'$u$'},'Interpreter','latex')
    
    subplot(1,3,3)
    plot(xb,pb,'-','Color',auxpar.linecolor)
    grid on
    set(gca,'GridLineStyle',':')
    title('$p$','Interpreter','latex')
    pause(0.001)

end

% The exact Riemann Solver is taken from here: https://www.mathworks.com/matlabcentral/fileexchange/48734-riemannexact-p1-rho1-u1-p4-rho4-u4-tol
Nex = 5*Nx;
[RHO,U,P] = RiemannExact(pL,rhoL,uL,pR,rhoR,uR,tend,xL,xR,Nex,1e-12);
X = linspace(xL,xR,Nex);

subplot(1,3,1)
hold on
plot(X,RHO)
hold off

subplot(1,3,2)
hold on
plot(X,U)
hold off

subplot(1,3,3)
hold on
plot(X,P)
hold off







