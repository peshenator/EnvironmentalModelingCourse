function Au = MatVecProd_u(u,Tb,Tcor) %,nu,nu,nu,nupy)

global dt dx dy Nx Ny

fm = zeros(Nx+1,Ny);
fp = zeros(Nx+1,Ny);

nub = visc(Tb); % viscosity at the barycenters
% x-direction 
fm(1     ,:) = 0;   % f_{-1}ghost cell (outside of the domain)
fm(2:Nx+1,:) =-nub(1:Nx,:).*( u(2:Nx+1,:) - u(1:Nx,:) )/dx;
fp(1:Nx  ,:) = fm(2:Nx+1,:);
fp(Nx+1  ,:) = 0;   % f_{Nx+1}ghost cell (outside of the domain)

Au = dt/dx*( fp - fm );

nuc = visc(Tcor);    % viscosity at the corners
% y-direction
fm(:,1)      =-nuc(:,   1).*( u(:,   1) - 0           )/(dy/2);
fm(:,2:Ny)   =-nuc(:,2:Ny).*( u(:,2:Ny) - u(:,1:Ny-1) )/dy;
fp(:,1:Ny-1) = fm(:,2:Ny);
fp(:,Ny)     =-nuc(:,  Ny).*( 0         - u(:,Ny    ) )/(dy/2);

Au = u + Au + dt/dy*( fp - fm );

end