%f BTCS scheme for the 2D Richards equation using the 
% nested Newton method of Casulli & Zanolli DOI:10.1137/100786320
clear;
% close all;
clc
global alpha thetas thetar n m Ks psic 
global dx dz dt K Kx Kz Nx Nz di Hx KL g 
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
Nx   = 100;   % horizontal 
Nz   = 20;    % vertical 
dx   = (xR-xL)/Nx; % x mesh spacing 
dz   = (zR-zL)/Nz; % y mesh spacing 
x    = linspace(xL+dx/2,xR-dx/2,Nx); 
z    = linspace(zL+dz/2,zR-dz/2,Nz); 
tend = 1e6;   % set the final simulation time 
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
    % Compute volume and hydraulic conductivity at the cell-centers
    theta(1:Nz,:) = Theta(psi(1:Nz,:));         % porous medium
    theta(Nz+1,:) = Height(psi(Nz+1,:),bath);   % free surface layer
    K(1:Nz,:) = Kfun(psi(1:Nz,:));              % porous medium
    K(Nz+1,:) = Kfun(psi(Nz+1,:));              % free surface layer

    thetan = theta; 
    if(nt==1)
        thetasum=sum(sum(theta))*dx*dz; 
    end
    % Compute depth of the free-surface layer at the cell-interfaces:
    H = theta(Nz+1,:);
    Hx(1   ) = 0;
    Hx(2:Nx) = max(0,max(H(2:Nx),H(1:Nx-1)));
    Hx(Nx+1) = 0;

    % Compute hydraulic conductivities at the cell-interfaces 
    % inside of the porous medium:
    Kx(1:Nz,1   ) = 0;
    Kx(1:Nz,2:Nx) = max(K(1:Nz,2:Nx),K(1:Nz,1:Nx-1));
    Kx(1:Nz,Nx+1) = 0;
    % artificial hydraulic conductivity in the free surface:
    Kx(Nz+1,:)    = g*dt*Hx.^2./(Hx + dt*gamma + 1e-14);

    Kz(1     ,1:Nx) = 0;
    Kz(2:Nz+1,1:Nx) = max(K(2:Nz+1,1:Nx),K(1:Nz,1:Nx));
    Kz(Nz+2  ,1:Nx) = 0;
    %Kz(Nz+1,:)=0e-13;

    
    % nonlinear convection of the free surface flow
    Fu = u; 
    
    % compute right hand side and the linear part of the system 
    [a,b,c,rhs] = LinearPartZ(theta,Fu,Hx,Kz);
    
    % Now, we start with the nested Newton method
    tol = 1e-7;   
    % initial guess
    psi(1:Nz,1:Nx) = min(psi(1:Nz,1:Nx),psic);
    for iNewton=1:100 % outer Newton iterations
        % The task of the outer iterations is ONLY to 
        % linearize one of the two nonlinear functions q1 or q2 
        di = zeros(Nz+1,Nx); % set the derivative of the nonlinear function to zero 
        Mpsi = matop2D(psi);  
        % compute the true nonlinear function f(psi)
        f(1:Nz,:) = Theta(psi(1:Nz,:)   )    + Mpsi(1:Nz,:)  - rhs(1:Nz,:); 
        f(Nz+1,:) = Height(psi(Nz+1,:),bath) + Mpsi(Nz+1 ,:) - rhs(Nz+1,:); 
        
        outres=sqrt(sum(sum(f.*f))); % outer residual
        disp(strcat(' Outer iteration = ', int2str(iNewton), ', outres = ', num2str(outres))); 
        if(outres<tol)
            break % tolerance has been reached 
        end
        psik = psi;          % save the value at the current outer iteration         
        %psi = max(psi,psic); % initial guess for the inner iterations 
        for inner = 1:100
            di   = zeros(Nz+1,Nx);  
            Mpsi = matop2D(psi);   % linear part of the system 
            
            fk(1:Nz,:) = Theta1(psi(1:Nz,:)) - (Theta2(psik(1:Nz,:)) + dTheta2(psik(1:Nz,:)).*(psi(1:Nz,:)-psik(1:Nz,:))) + Mpsi(1:Nz,:) - rhs(1:Nz,:);  
            di(1:Nz,:) = dTheta1(psi(1:Nz,:)) - dTheta2(psik(1:Nz,:)); 
            
            fk(Nz+1,:) = Height(psi(Nz+1,:),bath) + Mpsi(Nz+1,:) - rhs(Nz+1,:); % -(V2(psik(k,i))+dV2(psik(k,i))*(psi(k,i)-psik(k,i)))
            di(Nz+1,:) = dHeight(psi(Nz+1,:),bath);                      % -dTheta2(psik(k,i))

            for i=1:Nx
                for k=1:Nz+1
                    if(Kx(k,i) + Kx(k,i+1) + Kz(k,i) + Kz(k+1,i) == 0)
                        fk(k,i)=0;
                        di(k,i)=1; 
                    end
                end
            end

            inres = sqrt(sum(sum(fk.*fk))); 
            disp(strcat('  - Inner iteration = ',int2str(inner), ', inres = ', num2str(inres))) 
            if(inres<tol)
                break
            end
            [dpsi,CGres,CGk] = CGop2Dpc(fk);   % inner Newton step
            psi = psi - dpsi;      % update temperature at the inner iteration 
        end
    end 
    % now update the horizontal velocities
    u(2:Nx) = (Fu(2:Nx) - g*dt/dx*(psi(Nz+1,2:Nx)-psi(Nz+1,1:Nx-1)))./( 1+dt*gamma./Hx(2:Nx) );

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