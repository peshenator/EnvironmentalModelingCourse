% FTCS scheme for the nonlinear parabolic equation (Stefan problem)

clear;

global KS KL rhoS rhoL CS CL Tc TS TL hs kappaS kappaL;

% physical parameters
KS = 2.09;      % thermal conductivity of the ice (solid)
KL = 0.6;       % thermal conductivity of water
rhoS = 917;     % ice density
rhoL = 1000;    % water density
CS   = 2108;        % heat capacity of ice
CL   = 4187;        % heat capacity of water
kappaS   = KS/(rhoS*CS);
kappaL   = KL/(rhoL*CL);
Tc   = 0;           % CLitical temperature of the phase change
TS   = -10;         % left (air) temperature
TL   = +4;          % initial temperature of the lake
hs = 334e3;         % latent heat
day  = 24*3600;     % day in SI uniTc

% computational domain
xL = 0; 
xR = 2;
tend = 3*day;
Nx = 100;
dx = (xR - xL)/Nx;
x = linspace(xL+dx/2,xR-dx/2,Nx);
d = 0.45;       % d = K*dt/dx^2 < 0.5

% Initial conditions
T = TL*ones(1,Nx);
Q = zeros(1,Nx);
for i = 1:Nx
    Q(i) = Energy(T(i));
end

Km = zeros(1,Nx);
Kp = zeros(1,Nx);
fm = zeros(1,Nx);
fp = zeros(1,Nx);

% time loop
nmax = 100000;
time = 0;
for n = 1:nmax
    dt = d*dx^2/max(kappaS,kappaL);
    if (time + dt >tend)
        dt = tend - time;
    end
    if (time >= tend)
        break
    end
    
    % space loop
    Kb = K(T);
    % compute the thermal conductivity at the interfaces (Km at x_{i-1/2} and Kp at x_{i+1/2}})
    Km(1     ) = K(TS);
    Km(2:Nx  ) = 0.5*( Kb(2:Nx) + Kb(1:Nx-1) );
    Kp(1:Nx-1) = Km(2:Nx);
    Kp(Nx    ) = K(TL);

    fm(1     ) =-Km(1   )* ( T(1   ) - TS        )/(dx/2);
    fm(2:Nx  ) =-Km(2:Nx).*( T(2:Nx) - T(1:Nx-1) )/dx;
    fp(1:Nx-1) = fm(2:Nx);
    fp(Nx    ) =-Kp(Nx)*   ( TL      - T(Nx)     )/(dx/2);

    Q = Q - dt/dx*( fp - fm );
    
    for i=1:Nx
        T(i) = Temperature(Q(i));
    end
    
    time = time +dt;
    plot(x,T,'o-','Linewidth',1)
    title(strcat('Temperature at time = ',num2str(time)))
    grid on
    pause(0.01)
    
end

%% Compute the exact solution
hold on % add the exact solution to the existing plot 
% 
gamma = Newton(1);   % compute the constant gamma
% plot the exact solution at a given time
EMAX = 1000;                % number of poinTc to visualize the exact solution
xe = linspace(xL,xR,EMAX);  % set of poinTc where to visualize the exact solution
s = 2*gamma*sqrt(kappaS*time);  % location of the solid-liquid interface
G1 = (Tc-TS)/erf(gamma);    % temporary constanTc
G2 = (Tc-TL)/erfc(gamma*sqrt(kappaS/kappaL));
for i=1:EMAX
    if(xe(i)<s)
        % left branch (ice)
        Te(i) = TS + G1*erf( xe(i)/(2*sqrt(kappaS*time)) );
    else
        % right branch (liquid)
        Te(i) = TL + G2*erfc( xe(i)/(2*sqrt(kappaL*time)) );
    end
end
plot(xe,Te,'r-')
hold on
plot(s,Tc,'o','MarkerFaceColor',[1 0 0])
hold off
legend('FTCS','exact','interface')


















