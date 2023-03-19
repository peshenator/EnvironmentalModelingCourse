% Add the gradients of P' (Pprime) to the momemeta
function [ustar,vstar] = MomentumPressureGrad(ustar,vstar,Pstar)
global Nx Ny dx dy dt;

dtdx = dt/dx;
dtdy = dt/dy;

ustar(1   ,:) = ustar(1   ,:) - dtdx*(Pstar(1,:)      - Pstar(Nx,  :));
ustar(Nx+1,:) = ustar(Nx+1,:) - dtdx*(Pstar(1,:)      - Pstar(Nx,  :));
ustar(2:Nx,:) = ustar(2:Nx,:) - dtdx*(Pstar(2:Nx,:) - Pstar(1:Nx-1,:));


vstar(:,1   ) = vstar(:,1   ) - dtdy*(Pstar(:,1)      - Pstar(:,Ny  ));
vstar(:,Ny+1) = vstar(:,Ny+1) - dtdy*(Pstar(:,1)      - Pstar(:,Ny  ));
vstar(:,2:Ny) = vstar(:,2:Ny) - dtdy*(Pstar(:,2:Ny) - Pstar(:,1:Ny-1));

end
