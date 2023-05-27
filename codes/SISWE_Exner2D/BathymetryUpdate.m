% update a qb quantity (bathymetry in this case) at the barycenter using the Rusanov flux
function qb = BathymetryUpdate(qb,u,v,Hb)

global Nx Ny dt dx dy;

ub = 0.5*( u(2:Nx+1,:) + u(1:Nx,:) );
vb = 0.5*( v(:,2:Ny+1) + v(:,1:Ny) );

[fBx,fBy] = BathFlux(ub,vb,Hb);    % evaluate physical flux at the barecenters

% compute the numerical fluxes (Rusanov flux) at the cell faces
fx = zeros(Nx+1,Ny);
fy = zeros(Nx,Ny+1);
[aBx,aBy] = aBath(u,v,Hb);   % characterisitc velocity of the sediment at the cell faces
% Rusanov fluxes at the X and Y faces
fx(2:Nx,:) = 0.5*( fBx(2:Nx,:) + fBx(1:Nx-1,:) ) - 0.5*abs(aBx(2:Nx,:)).*(qb(2:Nx,:) - qb(1:Nx-1,:) );
fy(:,2:Ny) = 0.5*( fBy(:,2:Ny) + fBy(:,1:Ny-1) ) - 0.5*abs(aBy(:,2:Ny)).*(qb(:,2:Ny) - qb(:,1:Ny-1) );

% left X BC, periodic
fx(1   ,:) = 0.5*( fBx(1,:) + fBx(Nx,:) ) - 0.5*abs(aBx(1 ,:)).*( qb(1,:) - qb(Nx,:) );
% right X BC, periodic
fx(Nx+1,:) = 0.5*( fBx(1,:) + fBx(Nx,:) ) - 0.5*abs(aBx(Nx,:)).*( qb(1,:) - qb(Nx,:) );
% left  Y BC, wall (no flux)
fy(:,1   ) = fBy(:,1);
% right Y BC, wall (no flux)
fy(:,Ny+1) = fBy(:,Ny);


% finite volume update
qb = qb - dt/dx*( fx(2:Nx+1,:) - fx(1:Nx,:) ) - dt/dy*( fy(:,2:Ny+1) - fy(:,1:Ny) );

end