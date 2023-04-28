function p = pressure(rhox,vv,rhoE)
global gam Nx;

% kinetic energy at the integer points
rhoK = 0.25*( rhox(2:Nx+1).*vv(2:Nx+1).^2 + rhox(1:Nx).*vv(1:Nx).^2 );

p = (gam-1)*( rhoE - rhoK );

end
