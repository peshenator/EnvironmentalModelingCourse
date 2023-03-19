% Add the gradients of P' (Pprime) to the momemeta
function [u,v] = DivFreeVelocityCorrection(u,v,ustar,vstar,Pprime)
global Nx Ny dx dy dt;

dtdx = dt/dx;
dtdy = dt/dy;

u(1   ,:) = ustar(1   ,:) - dtdx*(Pprime(1   ,:) - Pprime(Nx,  :));
u(Nx+1,:) = ustar(Nx+1,:) - dtdx*(Pprime(1   ,:) - Pprime(Nx,  :));
u(2:Nx,:) = ustar(2:Nx,:) - dtdx*(Pprime(2:Nx,:) - Pprime(1:Nx-1,:));


v(:,1   ) = vstar(:,1   ) - dtdy*(Pprime(:,1   ) - Pprime(:,Ny  ));
v(:,Ny+1) = vstar(:,Ny+1) - dtdy*(Pprime(:,1   ) - Pprime(:,Ny  ));
v(:,2:Ny) = vstar(:,2:Ny) - dtdy*(Pprime(:,2:Ny) - Pprime(:,1:Ny-1));

end

