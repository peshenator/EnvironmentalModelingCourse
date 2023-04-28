function rhoE = energy(rhox,ux,p)
global gam Nx;

rho = 0.5*( rhox(2:Nx+1) + rhox(1:Nx));
u   = 0.5*( ux(2:Nx+1)   + ux(1:Nx));

% internal + kinetic
rhoE = p/(gam-1) + 0.5*rho.*u.^2;


