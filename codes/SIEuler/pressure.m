function p = pressure(rho,u,rhoE)
global gam;

% kinetic energy at the integer points
rhoK = 0.5*rho.*u.^2;

p = (gam-1)*( rhoE - rhoK );

end
