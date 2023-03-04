% Result of the matrix-vector product for the Poisson equation

function Ap=MultOpt2D_p(p)
global imax jmax dx dy

Ap = zeros(size(p));

% x-direction 
Ap(1,:)        = 1/dx*( (p(2,:)     -p(1,:)       )/dx     - (p(1,:)       -0            )/(dx/2) ); % Left BC
Ap(imax,:)     = 1/dx*( (0          -p(imax,:)    )/(dx/2) - (p(imax,:)    -p(imax-1,:)  )/dx ); % Right BC
Ap(2:imax-1,:) = 1/dx*( (p(3:imax,:)-p(2:imax-1,:))/dx     - (p(2:imax-1,:)-p(1:imax-2,:))/dx );

% y-direction 
Ap(:,1)    = Ap(:,1)    + 1/dy*( (p(:,2)-p(:,1))/dy - (p(:,1)-0)/(dy/2) ); % Bottom BC
Ap(:,jmax) = Ap(:,jmax) + 1/dy*( (0-p(:,jmax))/(dy/2) - (p(:,jmax)-p(:,jmax-1))/dy ); % Top BC
Ap(:,2:jmax-1) = Ap(:,2:jmax-1) + 1/dy*( (p(:,3:jmax)-p(:,2:jmax-1))/dy - (p(:,2:jmax-1)-p(:,1:jmax-2))/dy ); 

end
