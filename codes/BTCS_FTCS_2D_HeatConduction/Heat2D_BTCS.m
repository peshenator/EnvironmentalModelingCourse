%  BTCS for the 2D heat equation with the matrix-free Conjugate gradient method
% for solving the linear system
clear;
close all;
global dt dx dy Nx Ny;
global BCDirichletXL BCDirichletXR BCDirichletYL BCDirichletYR;

% numerical parameters
Nx = 64; Ny = Nx;
xL =-1; xR = 1; yL = -1; yR = 1;  % domain
t    = 0;
tend = 1;
d = 145; % for the FTCS scheme d <= 1/2
% BC parameters
BCDirichletXL = true; BCDirichletXR = true; BCDirichletYL = false; BCDirichletYR = false;
TxL = 80*ones(1,Ny);
TxR =  20*ones(1,Ny);
TyL = linspace(TxL(1),TxR(1),Nx);
TyR = linspace(TxL(1),TxR(1),Nx);

dx = (xR - xL)/Nx; dy = (yR - yL)/Ny;  % grid size
x = linspace(xL + dx/2, xR - dx/2, Nx);  % grid points x direction
y = linspace(yL + dy/2, yR - dy/2, Ny);  % grid points y direction

% physical parameters
L = zeros(Nx,Ny);
Lxm = zeros(Nx,Ny);
Lxp = zeros(Nx,Ny);
Lym = zeros(Nx,Ny);
Lyp = zeros(Nx,Ny);
[L,Lxm,Lxp,Lym,Lyp] = lambda(L,Lxm,Lxp,Lym,Lyp,x,y);  % thermal diffusivity lambda = K/(rho*C)
Lmax = max(max(L));

fm = zeros(Nx,Ny);
fp = zeros(Nx,Ny);
gm = zeros(Nx,Ny);
gp = zeros(Nx,Ny);


% initial condition
T = zeros(Nx, Ny);
for i=1:Nx
    for j=1:Ny
        if (x(i) < (x(1) + x(end))/2)
            T(i,j) = TxL(1);
        else
            T(i,j) = TxR(1);
        end
    end
end

%  plot initial condition
surf(x,y,T','EdgeColor','none','FaceColor','interp');
% view([0 90])

% time loop
for n = 1:10000000
    dt = d/(Lmax/dx^2 + Lmax/dy^2);
    if (t + dt > tend)
       dt = tend - t;
   end
   if (t >= tend)
       break;
   end

    % implementing BC in the BTCS scheme
    rhs = T;    % for Newuman BC rhs = T
    if (BCDirichletXL)
        rhs(1 ,:) = rhs(1 ,:) + dt/dx^2*Lxm(1,:).*TxL;
    end
    if (BCDirichletXR)
        rhs(Nx,:) = rhs(Nx,:) + dt/dx^2*Lxp(Nx,:).*TxR;
    end
    if (BCDirichletYL)
        rhs(:,1 ) = rhs(:,1 ) + dt/dy^2*Lym(:,1).*TyL;
    end
    if (BCDirichletYR)
        rhs(:,Ny) = rhs(:,Ny) + dt/dy^2*Lyp(:,Ny).*TyR;
    end
    
    % BTCS scheme
    T = CGsolver(rhs,@MatVectProd_T,Lxm,Lxp,Lym,Lyp);
    % update time
    t = t + dt;

    %  plot the solution
    surf(x,y,T','EdgeColor','none','FaceColor','interp');
    view([0 90])
    title(strcat('t = ',num2str(t),', n = ',num2str(n)))
    axis([xL xR yL yR -20 100])
    % clim([20 100])
    pause(0.001)
end


%% Compute the exact solution fo infinite domain (canbe used only for small times)
% problem (valid only for Riemann type initial data and not large times)
lambdaL = L(1 ,1);
lambdaR = L(Nx,1);
TL = TxL(1);
TR = TxR(1);
Tc = ( sqrt(lambdaR)*TR + sqrt(lambdaL)*TL )/( sqrt(lambdaL)+sqrt(lambdaR) );
imaxe=5*Nx;
xe = linspace(xL,xR,imaxe); % very fine grid for the plot of the exact solution 
Te = zeros(1,imaxe);
for i=1:imaxe
    if(xe(i)<0)
        Te(i) = TL+(Tc-TL)*erfc( -xe(i)/(2*sqrt(lambdaL*t)) ); 
    else
        Te(i) = TR+(Tc-TR)*erfc( +xe(i)/(2*sqrt(lambdaR*t)) ); 
    end
end
fig2 = figure;
plot(x,T(:,Ny/2),'o-')
hold on 
plot(xe,Te,'r-')
hold off
legend('FTCS','exact')

