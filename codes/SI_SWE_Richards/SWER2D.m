%f BTCS scheme for the 2D Richards equation using the 
% nested Newton method of Casulli & Zanolli DOI:10.1137/100786320
clear;
% close all;
clc
global alpha thetas thetar n m Ks psic 
global dx dz dt K Kx Kz Nx Nz di Hx KL g gamma
global a b c 
% Physical model parameters in SI units 
day     = 24*3600;    
g       = 9.81;         % m/s^2 
Ks      = 0.062/day;    % [m/s] 
thetas  = 0.41;         % [-] saturated water content 
thetar  = 0.095;        % [-] residual water content 
n       = 1.31;         % [-] parameter n 
m       = 1 - 1/n;      % [-] parameter m 
alpha   = 1.9;          % [m^(-1)] 
gamma   = 1e-2;         % friction coefficient 
% critical value of psi where the maximum of the C=dTheta/dPsi is located 
psic    = -1/alpha*( (n-1)/n )^(1/n); 
% Domain 
xL   = -1000; % left
xR   = +1000; % right 
zL   = -2;    % bottom 
zR   = 0;     % surface
Nx   = 200;   % horizontal 
Nz   = 20;    % vertical 
dx   = (xR-xL)/Nx; % x mesh spacing 
dz   = (zR-zL)/Nz; % y mesh spacing 
x    = linspace(xL+dx/2,xR-dx/2,Nx); 
z    = linspace(zL+dz/2,zR-dz/2,Nz); 
tend = 3000;   % set the final simulation time 
t    = 0;     % initial time 
% set the initial condition for pressure head
psi  = zeros(Nz+1,Nx);
psi(1:Nz,1:Nx) = -1;                                  % hydrostatic pressure
psi(Nz+1,1:Nx) = max(0, -0.5+exp(-0.5*x.^2/200^2) );  % free surface

u    = zeros(1,Nx+1);   % initial velocity of the free surface flow (at the cell faces) 
bath = zeros(1,Nx);     % bathymetry function (bottom shape)

K  = zeros(Nz+1,Nx);    % hydraulic conductivity at the cell centersin
Kx = zeros(Nz+1,Nx+1);  % hydraulic conductivity at the cell interfaces in x-direction
Kz = zeros(Nz+2,Nx);    % hydraulic conductivity at the cell interfaces in z-direction
Hx = zeros(1,Nx+1);     % free surface height at the cell interfaces
rhs= zeros(Nz+1,Nx);    % right hand side of the mildly nonlinear system

theta = zeros(Nz+1,Nx);
f     = zeros(Nz+1,Nx);    % true nonlinear function f(psi)
fk    = zeros(Nz+1,Nx);    
tic;
for nt=1:1000000000
    dt = 120/2;          % BTCS is unconditionally stable %0.45*dx^2/max(kL,kR); % time step restriction 
    if(t+dt>tend) 
        dt = tend-t; 
    end
    if(t>=tend)
        break
    end

    % nonlinear convection of the free surface flow
    Fu = u; 
    
    % compute right hand side and the linear part of the system 
    [a,b,c,theta,rhs] = LinearPartZ(psi,theta,bath,Fu);
    if(nt==1)
        thetasum=sum(sum(theta))*dx*dz; 
    end 
   
    % Now, we start with the nested Newton method
    tol = 1e-7;   
    % ------------ Outer Iterations -------------
    psi(1:Nz,1:Nx) = min(psi(1:Nz,1:Nx),psic);  % initial guess for the outer iterations
    for iOut=1:100 % outer Newton iterations
        % The task of the outer iterations is ONLY to 
        % linearize one of the two nonlinear functions Theta1 or Theta2 
        [ResOut,NormResOut] = ResidualOuter(psi,bath,rhs);
        disp(strcat(' Outer iteration = ', int2str(iOut), ', outres = ', num2str(NormResOut))); 
        if(NormResOut < tol)
            break
        end
        psik = psi;          % save the value at the current outer iteration         
        % ------ Inner Iterations ---------------        
        % psi = max(psi,psic); % initial guess for the inner iterations 
        for iInn = 1:100
            [ResInn,NormResInn] = ResidualInner(psi,psik,bath,rhs);
            disp(strcat('  - Inner iteration = ',int2str(iInn), ', inres = ', num2str(NormResInn))) 
            if(NormResInn < tol)
                break
            end
            [dpsi,CGres,CGk] = CGsolver(ResInn,@MatVecProd_psi);   % inner Newton step
            psi = psi - dpsi;
        end
    end 
    % now update the horizontal velocities
    u(2:Nx) = (Fu(2:Nx) - g*dt/dx*(psi(Nz+1,2:Nx) - psi(Nz+1,1:Nx-1)))./( 1 + dt*gamma./Hx(2:Nx) );

    % update time
    t = t + dt;

    % plotting and postprocessing
    masserror = thetasum - sum(sum(theta))*dx*dz 

    subplot(3,1,1)
    plot(x,theta(Nz+1,:),'o-')
    axis([xL xR 0 0.7])
    xlabel('x')
    ylabel('\eta')
    title(strcat('t = ',num2str(t))) 

    subplot(3,1,2) 
    surf(x,z,psi(1:Nz,:),'EdgeColor','none','FaceColor','interp')   
    xlabel('x')
    ylabel('z')
    zlabel('\psi')
    % view([0 90])
    title('\psi') 
    
    subplot(3,1,3) 
    surf(x,z,Theta(psi(1:Nz,:)),'EdgeColor','none','FaceColor','interp')   
    xlabel('x')
    ylabel('z')
    zlabel('\theta(\psi)')
    % view([0 90])
    title('\theta(\psi)') 
    
    pause(0.0001) 

end
toc