%  BTCS for the 2D heat equation with the matrix-free Conjugate gradient method
% for solving the linear system
clear;
global dt dx dy Nx Ny;

% numerical parameters
xL =-1; xR = 1; yL = -1; yR = 1;  % domain
t    = 0;
tend = 0.01;
d = 0.45; % for the FTCS scheme d <= 1/2

Nx = 64; Ny = Nx;  % number of grid points in x and y direction
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
            T(i,j) = 100;
        else
            T(i,j) = 20;
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

    % the FTCS scheme
    fm(1     ,:) = 0;  % boundary condition, no heat flux at the left boundary
    fm(2:Nx  ,:) =-Lxm(2:Nx,:).*( T(2:Nx,:) - T(1:Nx-1,:) )/dx;   % Fourier law
    fp(1:Nx-1,:) = fm(2:Nx,:);
    fp(Nx    ,:) = 0;  % boundary condition, no heat flux at the right boundary

    gm(:,1     ) = 0;  % boundary condition, no heat flux at the left boundary
    gm(:,2:Ny  ) =-Lym(:,2:Nx).*( T(:,2:Ny) - T(:,1:Ny-1) )/dy;   % Fourier law
    gp(:,1:Ny-1) = gm(:,2:Ny);
    gp(:,  Ny  ) = 0;  % boundary condition, no heat flux at the right boundary

    T = T - dt/dx*( fp - fm ) - dt/dy*( gp - gm );  % update temperature
    % update time
    t = t + dt;

    %  plot the solution
    surf(x,y,T','EdgeColor','none','FaceColor','interp');
    view([0 90])
    title(strcat('t = ',num2str(t),', n = ',num2str(n)))
    axis([xL xR yL yR 20 100])
    clim([20 100])
    pause(0.001)
end


%% Compute the exact solution fo infinite domain (canbe used only for small times)
% problem (valid only for Riemann type initial data and not large times)
lambdaL = L(1 ,1);
lambdaR = L(Nx,1);
TL = T(1,1);
TR = T(Nx,1);
Tc = ( sqrt(lambdaR)*TR + sqrt(lambdaL)*TL )/( sqrt(lambdaL)+sqrt(lambdaR) );
imaxe=5*Nx;
xe = linspace(xL + dx/2,xR - dx/2,imaxe); % very fine grid for the plot of the exact solution 
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

