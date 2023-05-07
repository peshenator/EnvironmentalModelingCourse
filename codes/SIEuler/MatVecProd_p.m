%  matrix-vector product for the pressure subsystem
function Mp = MatVecProd_p(p,hx,dt,dx)
global gam Nx;

fm = zeros(1,Nx);
fp = zeros(1,Nx);

dtdx = dt/dx;

fm(1     ) = hx(1   )*p(1 )*dtdx;   % hx(1).*(p(1) - pL)*dtdx but pL goes to the right hand side
fm(2:Nx  ) = hx(2:Nx).*(p(2:Nx)-p(1:Nx-1))*dtdx;
fp(1:Nx-1) = fm(2:Nx);
fp(Nx    ) =-hx(Nx+1)*p(Nx)*dtdx;% hx(Nx+1).*(pR - p(Nx))*dtdx but pR goes to the right hand side

Mp =-dtdx*(fp - fm) + p/(gam-1);

end