% Compute spatial differential operator of the energy equation
function [Km,Kp,rhs] = LinearPartCoeff(T)

global dt dx Nx TL TR;


Km = zeros(1,Nx);
Kp = zeros(1,Nx);

Km(1)      = K(TL);
Km(2:Nx)   = 0.5*( K(T(2:Nx)) + K(T(1:Nx-1)) );
Kp(1:Nx-1) = Km(2:Nx);
Kp(Nx)     = K(TR);

dtdx2 = dt/dx^2;

% compute the internal energy Q(T)
rhs(1     ) = Q(T(1     )) + dtdx2*2*Km(1)*TL;     
rhs(2:Nx-1) = Q(T(2:Nx-1));
rhs(Nx    ) = Q(T(Nx    )) + dtdx2*2*Kp(Nx)*TR;


end