close;
clear;
% BTCS scheme for the Stefan problem based on the nested Newton method of
% Casulli and Zanolli
global kappaS kappaL hs cS cL rhoL rhoS KS KL Tair Tlake Ts epsilon Nx Ny dt dx dy
% regularization parameter for the internal energy Q(T)
epsilon = 0.05;
% physical parameters of water and ice [in SI units]
KS = 2.09;     % heat conductivity of the solid phase (ice)
KL = 0.6;      % ... of the liquid phase (water)
hs = 334e3;         % specific latent heat
rhoS = 917;         % density of the solid
rhoL = 1000;        % density of the liquid
cS   = 2108;        % heat capacity of the solid
cL   = 4187;        % heat capacity of the liquid
kappaS   = KS/(rhoS*cS);
kappaL   = KL/(rhoL*cL);
Tc   =-0.1;           % critical temperature of the phase change
Ts   = Tc;
Tair = -10;         % top (air) temperature
Tlake   = 4;       % initial temperature of the lake
day  = 24*3600;     % day in SI units

% domain
xL = 0;             % left boundary
xR = 1;             % right boundary
yB = 0;             % bottom boundary (bottom of the lake)
yT = 1;             % top boundary (free surface)

tend = 3*day; %20*day;       % final time
% discretization parameters
Nx = 64;         % number of control volumes
Ny = Nx;         % number of control volumes
dx = (xR-xL)/Nx;  % mesh spacing in x
dy = (yT-yB)/Ny;  % mesh spacing in y

d  = 100;          % stability parameter for FTCS d = k/(dx^2/dt + dy^2/dt) < 0.5, or 0.25 in 2D ?
% mesh & initial condition
T = zeros(Nx,Ny);
T(:,:) = Tlake;
x = linspace(xL + dx/2,xR - dx/2,Nx);
y = linspace(yB + dy/2,yT - dy/2,Ny);
[X,Y] = meshgrid(x,y);

ResOut = zeros(Nx,Ny);             % Residual of the Outer iterations
ResIn  = zeros(Nx,Ny);             % Residual of the Inner iterations
dQ12   = zeros(Nx,Ny);             % diagonal elements of the Jacobian of Q1-Q2

time = 0;                           % initial time
nmax = 100000;                      % maximum number of time steps
tic
for n=1:nmax
    dt = d/( max(kappaS,kappaL)/dx^2 + max(kappaS,kappaL)/dy^2 );         % time step restriction of FTCS
    if(time+dt>tend)
        dt=tend-time;
    end
    if(time>=tend)
        break
    end
    % Compute the three diagonals of the matrix M
    [Kmx,Kpx,Kmy,Kpy,rhs] = LinearPartCoeff(T);
    T0 = T;
%%  Newsted Newton method of Casulli & Zanolli
    tol = 1e-12*rhoL*hs;
    T0 = min(T0,Ts-epsilon); % Initial guess for the outer iterations, see the paper by Casulli & Zanolli
    MaxNewton = 100;
    % ----- OUTER iterations -------
    for iouter = 1:MaxNewton
        ResOut = ResidualOuter(T0,Kmx,Kpx,Kmy,Kpy,rhs);    % Compute the residual of the outer iterations
        ResOut_norm = norm(ResOut);             % the norm of the residual
        if (ResOut_norm < tol)
            break
        end
        Talpha = T0; % we store the value of the temperature so that from now on the meaning of T0 is T^(alpha,beta)
        % ----- INNER iterations ------
        T0 = max(T0,Ts-epsilon);   % Initial guess for the inner iterations, see the paper by Casulli & Zanolli
        for iinner = 1:MaxNewton
            [ResInner,dQ12] = ResidualInner(T0,Talpha,Kmx,Kpx,Kmy,Kpy,rhs);    % Compute the residual of the inner iterations and the Jacobian of Q=Q1-Q2
            ResInner_norm = norm(ResInner);
%             disp(strcat('        Inner iteration = ',num2str(iinner),'|| residual =',num2str(ResInner_norm)))
            if (ResInner_norm < tol)
                break
            end
            % dT = Thomas(a,b+dQ12,c,ResInner);
            [dT,CGiter,CGerr] = CGsolver(ResInner,@MatVecProd_T,Kmx,Kpx,Kmy,Kpy,dQ12);
            T0 = T0 - dT;
        end
    end
    T = T0;
    time = time + dt; % advance time
    surf(x,y,T','FaceColor','interp','EdgeColor','none')     % plot the temperature profile 
    view([0 90])
    title(sprintf('Current time = %f',time))
    hold on
    contour(X,Y,T',[ 0 0 ],'b', 'LineWidth', 2)
    hold off
    pause(0.001)     
end
hold on % add the exact solution to the existing plot 
% 
gamma = Newton(1);   % compute the constant gamma
% plot the exact solution at a given time
EMAX = 1000;                % number of points to visualize the exact solution
Te = zeros(1,EMAX);         % exact solution
xe = linspace(xL,xR,EMAX);  % set of points where to visualize the exact solution
s = 2*gamma*sqrt(kappaS*time);  % location of the solid-liquid interface
G1 = (Tc-Tair)/erf(gamma);    % temporary constants
G2 = (Tc-Tlake)/erfc(gamma*sqrt(kappaS/kappaL));
for i=1:EMAX
    if(xe(i)<s)
        % left branch (ice)
        Te(i) = Tair + G1*erf( xe(i)/(2*sqrt(kappaS*time)) );
    else
        % right branch (liquid)
        Te(i) = Tlake + G2*erfc( xe(i)/(2*sqrt(kappaL*time)) );
    end
end
fig2 = figure;
plot(xe(end:-1:1),Te)
hold on
plot(y,T(Nx/2,:),'o-')
% plot(s,Ts,'o','MarkerFaceColor',[1 0 0])
legend('BTCS','exact')
hold off



toc







