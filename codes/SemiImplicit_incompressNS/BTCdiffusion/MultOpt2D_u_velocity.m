function Au = MultOpt2D_u_velocity(u)

global dt dx dy nu imax jmax uWall uLid;

Au = zeros(size(u));

% x-direction 
Au(1,:)        = - nu*dt/(dx)*( (u(2,:)     -u(1,:)       )/dx     - (u(1,:)       -0            )/(dx/2) );      % Left BC
Au(imax,:)     = - nu*dt/(dx)*( (0          -u(imax,:)    )/(dx/2) - (u(imax,:)    -u(imax-1,:)  )/dx );      % Right BC
Au(2:imax-1,:) = - nu*dt/(dx)*( (u(3:imax,:)-u(2:imax-1,:))/dx     - (u(2:imax-1,:)-u(1:imax-2,:))/dx );

% Au(imax+1,:) = - nu*dt/(dx)*( (0          -u(imax,:)    )/(dx/2) - (u(imax,:)    -u(imax-1,:)  )/dx );      % Right BC
% Au(2:imax,:) = - nu*dt/(dx)*( (u(3:imax+1,:)-u(2:imax,:))/dx     - (u(2:imax,:)-u(1:imax-1,:))/dx );


% y-direction 
Au(:,1)        = u(:,1) + Au(:,1)        - nu*dt/(dy)*( (u(:,2)     -u(:,1)       )/dy     - (u(:,1)       -0            )/(dy/2) );    % Bottom BC
Au(:,jmax)     = u(:,jmax) + Au(:,jmax)     - nu*dt/(dy)*( (0          -u(:,jmax)    )/(dy/2) - (u(:,jmax)    -u(:,jmax-1)  )/dy );    % Top BC
Au(:,2:jmax-1) = u(:,2:jmax-1) + Au(:,2:jmax-1) - nu*dt/(dy)*( (u(:,3:jmax)-u(:,2:jmax-1))/dy     - (u(:,2:jmax-1)-u(:,1:jmax-2))/dy );

 
% Au(:,1)        = Au(:,1)        - nu*dt/(dy)*( (u(:,2)     -u(:,1)       )/dy     - (u(:,1)       -0            )/(dy/2) );    % Bottom BC
% Au(:,jmax)     = Au(:,jmax)     - nu*dt/(dy)*( (0          -u(:,jmax)    )/(dy/2) - (u(:,jmax)    -u(:,jmax-1)  )/dy );    % Top BC
% Au(:,2:jmax-1) = Au(:,2:jmax-1) - nu*dt/(dy)*( (u(:,3:jmax)-u(:,2:jmax-1))/dy     - (u(:,2:jmax-1)-u(:,1:jmax-2))/dy );



end