function AT = MatVecProd_T(T)

global dt dx dy lambda Nx Ny;

fm = zeros(Nx,Ny);
fp = zeros(Nx,Ny);
gm = zeros(Nx,Ny);
gp = zeros(Nx,Ny);

% x direction fluxes
fm(1   ,:) = 0;                              % No heat flux (Neumann type BC)
fm(2:Nx  ,:) =- lambda*(T(2:Nx,:) - T(1:Nx-1,:))/dx;
fp(1:Nx-1,:) = fm(2:Nx  ,:);
fp(Nx    ,:) = 0; 

% y direction fluxes
gm(:,1     ) = 0; %0.5*v(:,1   ).*( T(:,1   ) + Tlake       ) - 0.5*abs(v(:,1   )).*( T(:,1   ) - Tlake       );
gm(:,2:Ny  ) =- lambda*(T(:,2:Ny) - T(:,1:Ny-1));
gp(:,1:Ny-1) = gm(:,2:Ny  );
gp(:,Ny    ) = 0;

AT = T + dt/dx*(fp - fm) + dt/dy*(gp - gm);  % Finite Volume update


end