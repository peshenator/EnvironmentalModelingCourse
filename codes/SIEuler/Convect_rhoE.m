
function rhoE_new = Convect_rhoE(rhoE,rhox,vx,ivar)
global Nx dt dx QL QR;

rhoE_new = rhoE;
dtdx = dt/dx;

% kinetic energy at the cell-centers:
rhoKc = 0.25*( rhox(2:Nx+1).*vx(2:Nx+1).^2 + rhox(1:Nx).*vx(1:Nx).^2 );

fx = zeros(1,Nx+1);
fx(2:Nx) = 0.5*( vx(2:Nx).*( rhoKc(2:Nx) + rhoKc(1:Nx-1) )- abs(vx(2:Nx)).*( rhoKc(2:Nx) - rhoKc(1:Nx-1) ) );

rhoE_new(1)      = QL(ivar);
rhoE_new(2:Nx-1) = rhoE(2:Nx-1) - dtdx*( fx(3:Nx) - fx(2:Nx-1) );
rhoE_new(Nx)     = QR(ivar);

end