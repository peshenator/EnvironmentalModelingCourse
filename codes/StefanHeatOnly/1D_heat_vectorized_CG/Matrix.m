% a is the lower diagonal
% b is the main diagonal
% c is the upper diagonal
% rhs is the right hand side of the mildly nonlinear system
function [a,b,c,rhs] = Matrix(T)

global dt dx Nx TL TR;

a   = zeros(1,Nx);                  % lower diagonal of the M matrix
b   = zeros(1,Nx);                  % main diagonal of the M matrix
c   = zeros(1,Nx);                  % upper diagonal of the M matrix
rhs = zeros(1,Nx);                  % upper diagonal of the M matrix
Kp  = zeros(1,Nx);
Km  = zeros(1,Nx);

dtdx2 = dt/dx^2;

% heat conductivity on the cell surfaces
Km(1)      = K(TL);
Km(2:Nx)   = 0.5*(K(T(2:Nx)) + K(T(1:Nx-1)));
Kp(1:Nx-1) = Km(2:Nx);
Kp(Nx)     = K(TR);

% elemnts of the tri-diagonal metriax
% inner points 
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

