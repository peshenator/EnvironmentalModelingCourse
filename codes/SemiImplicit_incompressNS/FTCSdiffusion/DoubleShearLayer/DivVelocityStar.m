% compute the right-hand side (rhs) for the pressure Poisson equation = 1/dt*div(u*,v*)
function divuv = DivVelocityStar(ustar,vstar)
global dx dy dt imax jmax;


ux = (ustar(2:imax+1,:) - ustar(1:imax,:))/dx;
vy = (vstar(:,2:jmax+1) - vstar(:,1:jmax))/dy;

divuv = (1/dt)*(ux + vy);

end