% rhs is the right hand side of the mildly nonlinear system, rhs = Q(T^n) + BC
function [Kmx,Kpx,Kmy,Kpy,rhs] = LinearPartCoeff(T)

global dt dx dy Nx Ny Tair Tlake;

rhs = zeros(Nx,Ny);                  % upper diagonal of the M matrix

dtdx2 = dt/dx^2;
dtdy2 = dt/dy^2;

% compute Km = K_{i-1/2} and Kp = K_{i+1/2} in each directions
Kmx = zeros(Nx,Ny);
Kpx = zeros(Nx,Ny);
Kmy = zeros(Nx,Ny);
Kpy = zeros(Nx,Ny);

Kmx(1     ,:) = K(Tlake);
Kmx(2:Nx  ,:) = 0.5*( K(T(2:Nx,:)) + K(T(1:Nx-1,:)) );
Kpx(1:Nx-1,:) = Kmx(2:Nx,:);
Kpx(Nx    ,:) = K(Tlake);

Kmy(:,1     ) = K(Tlake);
Kmy(:,2:Ny  ) = 0.5*( K(T(:,2:Ny)) + K(T(:,1:Ny-1)) );
Kpy(:,1:Ny-1) = Kmy(:,2:Ny);
Kpy(:,Ny    ) = K(Tair);


% inner points 
rhs(2:Nx-1,2:Ny-1) = Q(T(2:Nx-1,2:Ny-1));

% Left boundar
rhs(1,: ) = Q(T(1,: )) + dtdx2*2*Tlake*Kmx(1 ,:);
% Right BC
rhs(Nx,:) = Q(T(Nx,:)) + dtdx2*2*Tlake*Kpx(Nx,:);

% Bottom BC
rhs(:,1 ) = Q(T(:,1 )) + dtdy2*2*Tlake*Kmy(:,1 );
% Top BC (lake surface)
rhs(:,Ny) = Q(T(:,Ny)) + dtdy2*2*Tair *Kpy(:,Ny);

end

