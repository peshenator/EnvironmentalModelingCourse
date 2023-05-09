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

end