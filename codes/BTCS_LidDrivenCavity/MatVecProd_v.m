function Av = MatVecProd_v(v)

global dt dx dy nu Nx Ny vWall;

fm = zeros(Nx,Ny+1);
fp = zeros(Nx,Ny+1);

% y - direction
fm(:,1     ) =-nu*( v(:,1     ) - 0         )/dy;   % f_{-1}ghost cell (outside of the domain)
fm(:,2:Ny+1) =-nu*( v(:,2:Ny+1) - v(:,1:Ny) )/dy; 
fp(:,1:Ny  ) = fm(:,2:Ny+1);
fp(:,Ny+1  ) =-nu*( 0           - v(:,Ny+1) )/dy;   % f_{Nx+1}ghost cell (outside of the domain)

Av = dt/dy*( fp - fm );

% x - direction
fm(1     ,:) =-nu*( v(   1,:) - 0           )/(dx/2);
fm(2:Nx  ,:) =-nu*( v(2:Nx,:) - v(1:Nx-1,:) )/dx;
fp(1:Nx-1,:) = fm(2:Nx  ,:);
fp(Nx    ,:) =-nu*( 0         - v(Nx    ,:) )/(dx/2);

Av = v + Av + dt/dx*( fp - fm );

% Av = zeros(size(v));
% 
% % x-direction 
% Av(1       ,:) =-nu*dt/dx*( (v(2     ,:)-v(1       ,:))/dx     - (v(1       ,:) -vWall        )/(dx/2) );
% Av(2:Nx-1,:) =-nu*dt/dx*( (v(3:Nx,:)-v(2:Nx-1,:))/dx     - (v(2:Nx-1,:) -v(1:Nx-2,:))/dx     );
% Av(Nx    ,:) =-nu*dt/dx*( (vWall      -v(Nx    ,:))/(dx/2) - (v(Nx    ,:) -v(Nx-1,:)  )/dx     );
% 
% % y-direction 
% Av(:,1)      = v(:,1     ) + Av(:,1     ) - nu*dt/dy*( (v(:,2)-v(:,1) )/dy     - (v(:,1     )-0            )/(dy/2));
% Av(:,2:Ny) = v(:,2:Ny) + Av(:,2:Ny) - nu*dt/dy*( (v(:,3:Ny+1)-v(:,2:Ny))/dy     - (v(:,2:Ny)-v(:,1:Ny-1))/dy    );
% Av(:,Ny+1) = v(:,Ny+1) + Av(:,Ny+1) - nu*dt/dy*( ( (0-v(:,Ny+1)           )/(dy/2) - (v(:,Ny+1)-v(:,Ny)    )/dy )   );


end