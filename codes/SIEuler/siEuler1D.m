clear;
% close all;
global Nx dx dt cv rho0 gam uL uR pL pR QL QR plotcall xL xR ;
global MaxPicardP tolPicard numFlux MUSCL;

% material parameters for Ideal gas Equation of State
cv = 1;
gam = 1.4;
rho0 = 1;

[rhoL, uL, pL, rhoR, uR, pR, xL, xR, tend, Nd] = RPinit(2);

% discretization parameters
numFlux = 1;    % 1 = Rusanov, 2 = FORCE
MUSCL   = 0;    % 0 = noMUSCL, 1 = MUSCL
Nx = 1024;
dx = (xR - xL)/Nx;
xb = linspace(xL+dx/2,xR-dx/2,Nx  ); % cell- centers
x  = linspace(xL,xR,Nx+1); % cell edges
CFL = 0.5;

MaxPicardP = 50;
tolPicard  = 1e-6;

%% Initial data and BC: 
% cell x-face quantities:
qx = zeros(3,Nx+1);

qx(1,1:Nd) = rhoL;        qx(1,Nd:Nx+1) = rhoR;    
qx(2,1:Nd) = rhoL*uL;     qx(2,Nd:Nx+1) = rhoR*uR;
px(  1:Nd) = pL;          px(  Nd:Nx+1) = pR;
qx(3,:) =  px/(gam-1) + 0.5*qx(2,:).^2./qx(1,:); % rhoE_x(qx,px);

% barycenter quantities:
qb = 0.5*(qx(:,2:Nx+1) + qx(:,1:Nx));
pb = 0.5*(px(2:Nx+1)   + px(1:Nx));

QL = [rhoL rhoL*uL qx(3,1   )]';
QR = [rhoR rhoR*uR qx(  3,Nx+1)]';

time = 0;
plotcall = 0;
auxpar = struct('res',0,'linecolor','#4DBEEE');
myplot(time,xb,x,qx(1,:),qx(2,:)./qx(1,:),pb,qb(3,:),auxpar);

for n=1:10000000
    v_max =  max(abs(qx(2,:)./qx(1,:)));
    c_sound = sqrt(gam*pb./qb(1,:));
    c_max = max(v_max + c_sound);
    dtc  = CFL*dx/(c_max);
    if (n <= 50)
        dt = 3e-6;
    else
        dt  = CFL*dx/(v_max);
%         dt = 3e-3;
    end

    if (time >= tend)
        break
    end
    if (time + dt > tend)
        dt = tend  - time;
    end

    %% STEP #1 
    % explicit (convective sub-system, notation F(q) adapted from DOI:10.1016/j.amc.2015.08.042)
%     Fqb = Convect_Qb_MUSCL(qb,qx,dt,dx);   % convect all quantities at the cell-centeres
    Fqb = Convect_qb(qb,qx,dt,dx);   % convect all quantities at the cell-centeres
    Fqx = Convect_qx(qx,qb,dt,dx);   % convect all quantities at the cell-centeres
    
    %% STEP #2 
    % implicit (pressure sub-system)
    [qb,qx,residual,iter] = PressureSubSystem(Fqb,Fqx,dt,dx);
    % qb = Fqb; qx = Fqx;
     
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






