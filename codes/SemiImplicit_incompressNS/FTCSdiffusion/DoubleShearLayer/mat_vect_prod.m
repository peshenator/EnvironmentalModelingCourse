function Ap = mat_vect_prod(p)

global imax jmax dx dy
Ap = zeros(size(p));

% x-direction 
Ap(1,:)        = 1/dx^2*( (p(2,:)     -p(1,:)       ) - (p(1,:)       -p(imax,:)    ) ); % Left BC
Ap(imax,:)     = 1/dx^2*( (p(1,:)     -p(imax,:)    ) - (p(imax,:)    -p(imax-1,:)  ) ); % Right BC
Ap(2:imax-1,:) = 1/dx^2*( (p(3:imax,:)-p(2:imax-1,:)) - (p(2:imax-1,:)-p(1:imax-2,:)) );

% y-direction 
Ap(:,1)        = Ap(:,1)        + 1/dy^2*( (p(:,2)     -p(:,1)       ) - (p(:,1)        - p(:,jmax)    ) ); % Bottom BC
Ap(:,jmax)     = Ap(:,jmax)     + 1/dy^2*( (p(:,1)     -p(:,jmax)    ) - (p(:,jmax)     - p(:,jmax-1)  ) ); % Top BC
Ap(:,2:jmax-1) = Ap(:,2:jmax-1) + 1/dy^2*( (p(:,3:jmax)-p(:,2:jmax-1)) - (p(:,2:jmax-1) - p(:,1:jmax-2)) ); 

Ap = Ap + 1e-7*p;

end