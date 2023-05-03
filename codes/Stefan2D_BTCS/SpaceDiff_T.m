% Compute spatial differential operator of the energy equation
function DxyT = SpaceDiff_T(T,Kmx,Kpx,Kmy,Kpy)

global dt dx dy Nx Ny;

fm = zeros(Nx,Ny);
fp = zeros(Nx,Ny);

% X - direction
fm(1,:)      =-Kmx(1     ,:).*( T(1   ,:) - 0           )/(dx/2);
fm(2:Nx  ,:) =-Kmx(2:Nx  ,:).*( T(2:Nx,:) - T(1:Nx-1,:) )/dx;     % discrete Fourier law
fp(1:Nx-1,:) =-Kpx(1:Nx-1,:).*( T(2:Nx,:) - T(1:Nx-1,:) )/dx;     % discrete Fourier law
fp(Nx,:)     =-Kpx(Nx    ,:).*( 0         - T(Nx    ,:) )/(dx/2);   % No heat influx through the right boundary

DxyT = (dt/dx)*( fp - fm );

% Y - direction
fm(:,1)      =-Kmy(:,     1).*( T(:,   1) - 0           )/(dy/2);
fm(:,2:Ny)   =-Kmy(:,  2:Ny).*( T(:,2:Ny) - T(:,1:Ny-1) )/dy;     % discrete Fourier law
fp(:,1:Ny-1) =-Kpy(:,1:Ny-1).*( T(:,2:Ny) - T(:,1:Ny-1) )/dy;     % discrete Fourier law
fp(:,Ny)     =-Kpy(:,    Ny).*( 0         - T(:,Ny    ) )/(dy/2); % Dirichlet BC on the free surface

DxyT = DxyT + (dt/dy)*( fp - fm );

end