function [Lu,Lv] = ViscDiffOperator(u,v)

    global dt dx dy nu Nx Ny uWall vWall;
    
    Lu = zeros(size(u));
    Lv = zeros(size(v));
    
    % x-direction 
    Lu(1,:)        = - nu*dt/(dx)*( (u(2,:)     -u(1,:)       )/dx     - (u(1,:)       -uWall        )/(dx/2) );  % Left BC
    Lu(Nx,:)     = - nu*dt/(dx)*( (uWall      -u(Nx,:)    )/(dx/2) - (u(Nx,:)    -u(Nx-1,:)  )/dx );      % Right BC
    Lu(2:Nx-1,:) = - nu*dt/(dx)*( (u(3:Nx,:)-u(2:Nx-1,:))/dx     - (u(2:Nx-1,:)-u(1:Nx-2,:))/dx );
    
    % y-direction 
    Lu(:,1)        = Lu(:,1)        - nu*dt/(dy)*( (u(:,2)     -u(:,1)       )/dy     - (u(:,1)       -0            )/(dy/2) );% Bottom BC
    Lu(:,Ny)     = Lu(:,Ny)     - nu*dt/(dy)*( (0          -u(:,Ny)    )/(dy/2) - (u(:,Ny)    -u(:,Ny-1)  )/dy );    % Top BC
    Lu(:,2:Ny-1) = Lu(:,2:Ny-1) - nu*dt/(dy)*( (u(:,3:Ny)-u(:,2:Ny-1))/dy     - (u(:,2:Ny-1)-u(:,1:Ny-2))/dy );

    % x-direction 
    Lv(1,:)        = -nu*dt/(dx)*( (v(2,:)     -v(1,:)       )/dx     - (v(1,:)       -vWall        )/(dx/2) );      % Left BC
    Lv(Nx,:)     = -nu*dt/(dx)*( (vWall      -v(Nx,:)    )/(dx/2) - (v(Nx,:)    -v(Nx-1,:)  )/dx );          % Right BC
    Lv(2:Nx-1,:) = -nu*dt/(dx)*( (v(3:Nx,:)-v(2:Nx-1,:))/dx     - (v(2:Nx-1,:)-v(1:Nx-2,:))/dx );
    
    % y-direction 
    Lv(:,1)        = Lv(:,1)        - nu*dt/(dy)*( (v(:,2)-v(:,1)            )/dy     - (v(:,1)-0)/(dy/2));           % Bottom BC
    Lv(:,Ny)     = Lv(:,Ny)     - nu*dt/(dy)*( ((0-v(:,Ny)            )/(dy/2) - (v(:,Ny)-v(:,Ny-1))/dy)) ; % Top BC
    Lv(:,2:Ny-1) = Lv(:,2:Ny-1) - nu*dt/(dy)*( (v(:,3:Ny)-v(:,2:Ny-1))/dy     - (v(:,2:Ny-1)-v(:,1:Ny-2))/dy );
      
end