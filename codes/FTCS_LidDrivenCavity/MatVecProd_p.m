% Result of the matrix-vector product for the Poisson equation

function Ap=MatVecProd_p(p)
global Nx Ny dx dy

Ap = zeros(size(p));

% x-direction 
Ap(1     ,:) = 1/dx*( (p(2   ,:) - p(1     ,:))/dx     - (p(1,     :) - 0          )/(dx/2) ); % Left BC
Ap(Nx    ,:) = 1/dx*( (0         - p(Nx    ,:))/(dx/2) - (p(Nx,    :) - p(Nx-1  ,:))/dx ); % Right BC
Ap(2:Nx-1,:) = 1/dx*( (p(3:Nx,:) - p(2:Nx-1,:))/dx     - (p(2:Nx-1,:) - p(1:Nx-2,:))/dx );

% y-direction 
Ap(:,1     ) = Ap(:,1     ) + 1/dy*( (p(:,2   ) - p(:,1     ))/dy     - (p(:,1     ) - 0          )/(dy/2) ); % Bottom BC
Ap(:,Ny    ) = Ap(:,Ny    ) + 1/dy*( (0         - p(:,Ny    ))/(dy/2) - (p(:,Ny    ) - p(:,Ny-1)  )/dy ); % Top BC
Ap(:,2:Ny-1) = Ap(:,2:Ny-1) + 1/dy*( (p(:,3:Ny) - p(:,2:Ny-1))/dy     - (p(:,2:Ny-1) - p(:,1:Ny-2))/dy ); 

end
