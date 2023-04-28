% subroutine for the pressure sub-system
function [rhoE,mx,res,iPicard] = PSubsystemCG(rhox,mx,rhoE)

global Nx dt dx gamma pL pR QL QR;
maxPicard = 100;
tol = 1e-3;

hx = zeros(1,Nx+1);       % enthalpy at the x-faces
px = zeros(1,Nx+1);       % pressure at the x-faces

mstar = mx;                 % momentum after the advection step

% intial guess for the pressure
ux = mx./rhox;
p = pressure(rhox,ux,rhoE);    % pressure in the cell-centers

dtdx  = dt/dx;
dtdx2 = (dt/dx)^2;
res = 1;

% Picard loop (see the red "r" index in the lecture notes)
for iPicard = 1:maxPicard
    if (res<tol)
        break;
    end
    % update the coefficients of the matrix M (see the lecture notes)
    % update the enthealpy
    px(2:Nx) = 0.5*(p(2:Nx) + p(1:Nx-1)); % pressure at the vertexes
    px(1)    = pL;
    px(Nx+1) = pR;
    hx = px.*((gamma-1)*rhox) + px./rhox;
        
    % compute the right hand-side of the linear system
    ux   = mx./rhox;
    rhoK = 0.25*( rhox(2:Nx+1).*ux(2:Nx+1).^2 + rhox(1:Nx).*ux(1:Nx).^2 );
    rhs  = rhoE - rhoK  - dtdx*( hx(2:Nx+1).*mstar(2:Nx+1) - hx(1:Nx).*mstar(1:Nx) );
    rhs(1 ) = rhs(1 ) + dtdx2*hx(1   )*pL;
    rhs(Nx) = rhs(Nx) + dtdx2*hx(Nx+1)*pR;
    
    % solve the linear system
    res = norm(MatVecProdCG(p,hx,dt,dx) - rhs);
    p   = CGSolver(rhs,hx,dt,dx);
    
    % update the momentum
    mx(1)    = QL(2);
    mx(2:Nx) = mstar(2:Nx) - dtdx*( p(2:Nx) - p(1:Nx-1) );
    mx(Nx+1) = QR(2);
end

% update the total energy
rhoE = rhoE - dtdx*( hx(2:Nx+1).*mx(2:Nx+1) - hx(1:Nx).*mx(1:Nx) );



end