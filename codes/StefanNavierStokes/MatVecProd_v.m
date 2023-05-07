function Av = MatVecProd_v(v,Tb,Tcor)

global dt dx dy Nx Ny;

fm = zeros(Nx,Ny+1);
fp = zeros(Nx,Ny+1);

nub = visc(Tb);
% y - direction
fm(:,1     ) = 0;   % f_{-1}ghost cell (outside of the domain)
fm(:,2:Ny+1) =-nub(:,1:Ny).*( v(:,2:Ny+1) - v(:,1:Ny) )/dy; 
fp(:,1:Ny  ) = fm(:,2:Ny+1);
fp(:,Ny+1  ) = 0;   % f_{Nx+1}ghost cell (outside of the domain)

Av = dt/dy*( fp - fm );

nuc = visc(Tcor);
% x - direction
fm(1     ,:) =-nuc(   1,:).*( v(   1,:) - 0           )/(dx/2);
fm(2:Nx  ,:) =-nuc(2:Nx,:).*( v(2:Nx,:) - v(1:Nx-1,:) )/dx;
fp(1:Nx-1,:) = fm(2:Nx  ,:);
fp(Nx    ,:) =-nuc(  Nx,:).*( 0         - v(Nx    ,:) )/(dx/2);

Av = v + Av + dt/dx*( fp - fm );

end