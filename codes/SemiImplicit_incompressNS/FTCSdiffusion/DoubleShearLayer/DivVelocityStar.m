% compute the right-hand side (rhs) for the pressure Poisson equation = 1/dt*div(u*,v*)
function divuv = DivVelocityStar(ustar,vstar)
global dx dy dt Nx Ny;


ux = (ustar(2:Nx+1,:) - ustar(1:Nx,:))/dx;
vy = (vstar(:,2:Ny+1) - vstar(:,1:Ny))/dy;

divuv = (1/dt)*(ux + vy);

end