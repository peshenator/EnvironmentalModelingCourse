% Compute spatial differential operator of the energy equation
function DxyT = SpaceDiff_T(T,Km,Kp)

global dt dx Nx Ny TL TR;

fm = zeros(Nx,Ny);
fp = zeros(Nx,Ny);

% X - direction

fm(1)      =-Km(1   ).*( T(1) - 0 )/(dx/2);    % No heat influx through the left boundary
fm(2:Nx  ) =-Km(2:Nx  ).*( T(2:Nx) - T(1:Nx-1) )/dx;     % discrete Fourier law
fp(1:Nx-1) =-Kp(1:Nx-1).*( T(2:Nx) - T(1:Nx-1) )/dx;     % discrete Fourier law
fp(Nx)     =-Kp(Nx  ).*( 0 - T(Nx))/(dx/2);   % No heat influx through the right boundary

DxyT = (dt/dx)*( fp - fm );

% *******************************
a   = zeros(1,Nx);                  % lower diagonal of the M matrix
b   = zeros(1,Nx);                  % main diagonal of the M matrix
c   = zeros(1,Nx);                  % upper diagonal of the M matrix
rhs = zeros(1,Nx);                  % upper diagonal of the M matrix

% elemnts of the tri-diagonal metriax
% inner points 
dtdx2 = dt/dx^2;

a(2:Nx-1) =-dtdx2*  Km(2:Nx-1);
b(2:Nx-1) =+dtdx2*( Km(2:Nx-1) + Kp(2:Nx-1) );
c(2:Nx-1) =-dtdx2*  Kp(2:Nx-1);
rhs(2:Nx-1) = Q(T(2:Nx-1));        % compute the internal energy Q(T)

% Left boundar
a(1) = 0;
b(1) =+dtdx2*(2*Km(1) + Kp(1));
c(1) =-dtdx2*Kp(1);
rhs(1) = Q(T(1)) + dtdx2*2*Km(1)*TL;    % compute the internal energy Q(T)

% right BC
a(Nx) =-dtdx2*Km(Nx) ;
b(Nx) =+dtdx2*(Km(Nx) + 2*Kp(Nx));
c(Nx) = 0;
rhs(Nx) = Q(T(Nx)) + dtdx2*2*Kp(Nx)*TR;        % compute the internal energy Q(T)

end