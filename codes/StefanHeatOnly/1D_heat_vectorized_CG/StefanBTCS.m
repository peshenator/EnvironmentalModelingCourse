clear;
% BTCS scheme for the Stefan problem based on the nested Newton method of
% Casulli and Zanolli
global kappaL kappaR hs cL cR rhoL rhoR Ks Kl TL TR Ts epsilon Nx dt dx
% regularization parameter for the internal energy Q(T)
epsilon = 0.1;
% physical parameters of water and ice [in SI units]
Ks = 2.09;     % heat conductivity of ice
Kl = 0.6;      % ... of water
hs = 334e3;         % specific latent heat
rhoL = 917;         % density of ice
rhoR = 1000;        % density of water
cL   = 2108;        % heat capacity of ice
cR   = 4187;        % heat capacity of water
kappaL   = Ks/(rhoL*cL);
kappaR   = Kl/(rhoR*cR);
Tc   = 0;           % critical temperature of the phase change
Ts   = Tc;
TL   = -10;         % left (air) temperature
TR   = +4;          % initial temperature of the lake
day  = 24*3600;     % day in SI units
% domain
xL = 0;             % left boundary (free surface)
xR = 2;             % right boundary (bottom of the lake)
tend = 10*day; %20*day;       % final time
% discretization parameters
Nx = 64;         % number of control volumes
dx = (xR-xL)/Nx;  % mesh spacing
d  = 10;          % stability parameter for FTCS d = k*dt/dx^2 < 0.5
% mesh & initial condition
T = zeros(1,Nx);
for i=1:Nx
    x(i) = xL + dx/2 + (i-1)*dx;    % barycenters of the finite volumes
    T(i) = TR;                      % initial temperature
end

ResOut = zeros(1,Nx);             % Residual of the Outer iterations
ResIn  = zeros(1,Nx);             % Residual of the Inner iterations
dQ12   = zeros(1,Nx);             % diagonal elements of the Jacobian of Q1-Q2

time = 0;                           % initial time
nmax = 100000;                      % maximum number of time steps
tic
for n=1:nmax
    dt = d*dx^2/max(kappaL,kappaR);         % time step restriction of FTCS
    if(time+dt>tend)
        dt=tend-time;
    end
    if(time>=tend)
        break
    end
    % Construct the three diagonals of the matrix M
    [Km,Kp,rhs] = LinearPartCoeff(T);
    T0 = T;
%%  Newsted Newton method of Casulli & Zanolli
    tol = 1e-12*rhoR*hs;
    T0 = min(T0,Ts-epsilon);      % Initial guess for the outer iterations, see the paper by Casulli & Zanolli
    MaxNewton = 100;
    % ----- OUTER iterations -------
    for iouter = 1:MaxNewton
        ResOut = ResidualOuter(T0,Km,Kp,rhs);    % Compute the residual of the outer iterations
        ResOut_norm = norm(ResOut);             % the norm of the residual
%         disp(strcat('Outer iteration = ',num2str(iouter),'|| residual =',num2str(ResOut_norm)))
        if (ResOut_norm < tol)
            break
        end
        Talpha = T0; % we store the value of the temperature so that from now on the meaning of T is T^(alpha,beta)
        % ----- INNER iterations ------
        T0 = max(T0,Ts-epsilon);   % Initial guess for the inner iterations, see the paper by Casulli & Zanolli
        for iinner = 1:MaxNewton
            [ResInner,dQ12] = ResidualInner(T0,Talpha,Km,Kp,rhs);    % Compute the residual of the inner iterations and the Jacobian of Q=Q1-Q2
            ResInner_norm = norm(ResInner);
%             disp(strcat('        Inner iteration = ',num2str(iinner),'|| residual =',num2str(ResInner_norm)))
            if (ResInner_norm < tol)
                break
            end
            % dT = Thomas(a,b+dQ12,c,ResInner);
            [dT,CGiter,CGerr] = CG_T(ResInner,Km,Kp,dQ12);
            T0 = T0 - dT;
        end
    end
    T = T0;
    time = time + dt; % advance time
    plot(x,T,'-o','Linewidth',1)     % plot the temperature profile 
    title(sprintf('Current time = %f',time))
    grid on
    pause(0.001)     
end
hold on % add the exact solution to the existing plot 
% 
gamma = Newton(1)   % compute the constant gamma
% plot the exact solution at a given time
EMAX = 1000;                % number of points to visualize the exact solution
xe = linspace(xL,xR,EMAX);  % set of points where to visualize the exact solution
s = 2*gamma*sqrt(kappaL*time);  % location of the solid-liquid interface
G1 = (Tc-TL)/erf(gamma);    % temporary constants
G2 = (Tc-TR)/erfc(gamma*sqrt(kappaL/kappaR));
for i=1:EMAX
    if(xe(i)<s)
        % left branch (ice)
        Te(i) = TL + G1*erf( xe(i)/(2*sqrt(kappaL*time)) );
    else
        % right branch (liquid)
        Te(i) = TR + G2*erfc( xe(i)/(2*sqrt(kappaR*time)) );
    end
end
plot(xe,Te,'r-')
hold on
plot(s,Ts,'o','MarkerFaceColor',[1 0 0])
legend('BTCS','exact','ice-water interface')
hold off

toc







