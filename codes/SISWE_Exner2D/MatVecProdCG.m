function Teta = MatVecProdCG(etab,bathb,Hmx,Hpx,Hmy,Hpy,wet)

% This function computes the matrix-vector product T*eta
global g dt dx dy Nx Ny

kx = g*(dt/dx)^2;
ky = g*(dt/dy)^2;

fm = zeros(Nx,Ny);
fp = zeros(Nx,Ny);

% X - direction
fm(1     ,:) = Hmx(1     ,:).*( etab(1   ,:) - etab(  Nx  ,:) );
fm(2:Nx  ,:) = Hmx(2:Nx  ,:).*( etab(2:Nx,:) - etab(1:Nx-1,:) );
fp(1:Nx-1,:) = Hpx(1:Nx-1,:).*( etab(2:Nx,:) - etab(1:Nx-1,:) );
fp(  Nx  ,:) = Hpx(  Nx  ,:).*( etab(1   ,:) - etab(  Nx ,:)  );

Teta =-kx*( fp - fm );

% Y - direction
fm(:,1     ) = 0; %Hmy(:,     1).*( etab(:,   1) - 0              );
fm(:,2:Ny  ) = Hmy(:,2:Ny  ).*( etab(:,2:Ny) - etab(:,1:Ny-1) );
fp(:,1:Ny-1) = Hpy(:,1:Ny-1).*( etab(:,2:Ny) - etab(:,1:Ny-1) );
fp(:,  Ny  ) = 0; %Hpy(:,  Ny  ).*( 0            - etab(:,  Ny  ) );

Teta = Teta - ky*( fp - fm );

Teta = Teta + wet.*etab;

end