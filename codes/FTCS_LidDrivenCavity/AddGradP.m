% Add the gradients of P* (Pstar) to the momemeta
function [u,v] = AddGradP(u,v,u0,v0,P)
global Nx Ny dx dy dt;

dtdx = dt/dx;
dtdy = dt/dy;

% velocity correction v^{n+1} = v* - dt*grad(P'):
u(1   ,:) = 0;
u(2:Nx,:) = u0(2:Nx,:) - dtdx*(P(2:Nx,:) - P(1:Nx-1,:));
u(Nx+1,:) = 0;

v(:,1   ) = 0;
v(:,2:Ny) = v0(:,2:Ny) - dtdy*(P(:,2:Ny) - P(:,1:Ny-1));
v(:,Ny+1) = 0;


end
