% compute the result of the product (M + dTheta1 - dTheta2)*psi = M*psi + (dTheta1 - dTheta2)*psi
function Mpsi = MatVecProd_psi(psi)
global dx dz dt Kx Kz Nx Nz dTheta12
% diagonal part including the derivatives of the nonlinear functions 
% Mpsi = dTheta12.*psi;

fzm = zeros(Nz+1,Nx);
fzp = zeros(Nz+1,Nx);
fxm = zeros(Nz+1,Nx);
fxp = zeros(Nz+1,Nx);

% z fluxes
fzm(1     ,:) = 0;
fzm(2:Nz+1,:) = Kz(2:Nz+1 ,:).*( psi(2:Nz+1,:) - psi(1:Nz,:) )/dz;
fzp(1:Nz  ,:) = fzm(2:Nz+1,:);
fzp(Nz+1  ,:) = 0;
% <-- Explanation: for the free surface layer, there is no mass flux across the
% free surface, but a mass flux across the bottom, given by Richards equation ! 

%  x flux
fxm(:,1     ) = 0;
fxm(:,2:Nx  ) = Kx(:,2:Nx).*( psi(:,2:Nx) - psi(:,1:Nx-1) )/dx;
fxp(:,1:Nx-1) = fxm(:,2:Nx);
fxp(:,  Nx  ) = 0;

Mpsi =-dt/dz*( fzp - fzm ) - dt/dx*( fxp - fxm ) + dTheta12.*psi;

end