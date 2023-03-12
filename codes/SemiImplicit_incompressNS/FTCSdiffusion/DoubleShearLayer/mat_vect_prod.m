function Ap = mat_vect_prod(p)

global Nx Ny dx dy
Ap = zeros(size(p));

% x-direction 
Ap(1,:)        = 1/dx^2*( (p(2,:)     -p(1,:)       ) - (p(1,:)       -p(Nx,:)    ) ); % Left BC
Ap(Nx,:)     = 1/dx^2*( (p(1,:)     -p(Nx,:)    ) - (p(Nx,:)    -p(Nx-1,:)  ) ); % Right BC
Ap(2:Nx-1,:) = 1/dx^2*( (p(3:Nx,:)-p(2:Nx-1,:)) - (p(2:Nx-1,:)-p(1:Nx-2,:)) );

% y-direction 
Ap(:,1)        = Ap(:,1)        + 1/dy^2*( (p(:,2)     -p(:,1)       ) - (p(:,1)        - p(:,Ny)    ) ); % Bottom BC
Ap(:,Ny)     = Ap(:,Ny)     + 1/dy^2*( (p(:,1)     -p(:,Ny)    ) - (p(:,Ny)     - p(:,Ny-1)  ) ); % Top BC
Ap(:,2:Ny-1) = Ap(:,2:Ny-1) + 1/dy^2*( (p(:,3:Ny)-p(:,2:Ny-1)) - (p(:,2:Ny-1) - p(:,1:Ny-2)) ); 

Ap = Ap + 1e-7*p;

end