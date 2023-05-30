% compute the residual of the inner iterations and dTheta12 = dTheat1^{beta+1} - dTheta2^alpha
function [r,r_norm] = ResidualInner(psi,psik,bath,rhs)
global Nx Nz dTheta12

    r        = zeros(Nz+1,Nx);     % residual of the linearized system
    dTheta12 = zeros(Nz+1,Nx);  
    Mpsi     = MatVecProd_psi(psi);   % linear part of the system 
    
    r( 1:Nz,:) = Theta1(psi(1:Nz,:)) - Theta2(psik(1:Nz,:)) - dTheta2(psik(1:Nz,:)).*(psi(1:Nz,:)-psik(1:Nz,:)) + Mpsi(1:Nz,:) - rhs(1:Nz,:);  
    r( Nz+1,:) = Height(psi(Nz+1,:),bath) + Mpsi(Nz+1,:) - rhs(Nz+1,:);

    r_norm = sqrt(sum(sum(r.*r))); 

    dTheta12(1:Nz,:) = dTheta1(psi(1:Nz,:)) - dTheta2(psik(1:Nz,:)); 
    dTheta12(Nz+1,:) = dHeight(psi(Nz+1,:),bath);

end