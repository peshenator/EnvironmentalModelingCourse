function Av = MultOpt2D_v_velocity(v)

global dt dx dy nu imax jmax vWall;

Av = zeros(size(v));

% x-direction 
Av(1,:)        = -nu*dt/(dx)*( (v(2,:)     -v(1,:)       )/dx     - (v(1,:)       -vWall        )/(dx/2) );      % Left BC
Av(imax,:)     = -nu*dt/(dx)*( (vWall      -v(imax,:)    )/(dx/2) - (v(imax,:)    -v(imax-1,:)  )/dx );          % Right BC
Av(2:imax-1,:) = -nu*dt/(dx)*( (v(3:imax,:)-v(2:imax-1,:))/dx     - (v(2:imax-1,:)-v(1:imax-2,:))/dx );


% y-direction 
Av(:,1)        = v(:,1)        + Av(:,1)        - nu*dt/(dy)*( (v(:,2)-v(:,1)            )/dy     - (v(:,1)-0)/(dy/2));           % Bottom BC
Av(:,jmax)     = v(:,jmax)     + Av(:,jmax)     - nu*dt/(dy)*( ((0-v(:,jmax)            )/(dy/2) - (v(:,jmax)-v(:,jmax-1))/dy)) ; % Top BC
Av(:,2:jmax-1) = v(:,2:jmax-1) + Av(:,2:jmax-1) - nu*dt/(dy)*( (v(:,3:jmax)-v(:,2:jmax-1))/dy     - (v(:,2:jmax-1)-v(:,1:jmax-2))/dy );

end