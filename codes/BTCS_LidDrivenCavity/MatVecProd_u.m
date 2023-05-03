function Au = MatVecProd_u(u)

global dt dx dy nu Nx Ny uWall

fm = zeros(Nx+1,Ny);
fp = zeros(Nx+1,Ny);
uGhost = uWall;

% x-direction 
fm(1     ,:) =-nu*( u(1    ,:)  - uGhost    )/dx;      % f_{-1}ghost cell (outside of the domain)
fm(2:Nx+1,:) =-nu*( u(2:Nx+1,:) - u(1:Nx,:) )/dx;
fp(1:Nx  ,:) = fm(2:Nx+1,:);
fp(Nx+1  ,:) =-nu*( uGhost      - u(Nx+1,:) )/dx;       % f_{Nx+1}ghost cell (outside of the domain)

Au = dt/dx*( fp - fm );

% y-direction
fm(:,1)      =-nu*( u(:,   1) - 0           )/(dy/2);
fm(:,2:Ny)   =-nu*( u(:,2:Ny) - u(:,1:Ny-1) )/dy;
fp(:,1:Ny-1) = fm(:,2:Ny);
fp(:,Ny)     =-nu*( 0         - u(:,Ny    ) )/(dy/2);

Au = u + Au + dt/dy*( fp - fm );



% Au = zeros(size(u));
% 
% % % x-direction 
% Au(1,:)    =-nu*dt/(dx/2)*( (u(2,:) - u(1,:) )/dx - uWall );  % Left BC
% Au(2:Nx,:) =-nu*dt/dx    *( (u(3:Nx+1,:)-u(2:Nx,:))/dx - (u(2:Nx,:)-u(1:Nx-1,:))/dx );
% Au(Nx+1,:) =-nu*dt/(dx/2)*( uWall                          - (u(Nx+1,:)-u(Nx    ,:))/dx );  % Right BC
% 
% % y-direction 
% Au(:,1       ) = u(:,1       ) + Au(:,1       ) - nu*dt/dy*( (u(:,2)     -u(:,1       ))/dy     - (u(:,1       )-0        )/(dy/2) );% Bottom BC
% Au(:,2:Ny-1) = u(:,2:Ny-1) + Au(:,2:Ny-1) - nu*dt/dy*( (u(:,3:Ny)-u(:,2:Ny-1))/dy     - (u(:,2:Ny-1)-u(:,1:Ny-2))/dy     );
% Au(:,Ny    ) = u(:,Ny    ) + Au(:,Ny    ) - nu*dt/dy*( (0       -u(:,Ny    ))/(dy/2) - (u(:,Ny    )-u(:,Ny-1  ))/dy     );    % Top BC

end