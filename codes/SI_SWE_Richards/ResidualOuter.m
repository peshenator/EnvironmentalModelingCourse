% compute the residual of the outer iterations of the newsted Newton method
function [r,r_norm] = ResidualOuter(psi,bath,rhs)
global dTheta12 Nx Nz

    r        = zeros(Nz+1,Nx);      % residual of the original mildly nonlinear system
    dTheta12 = zeros(Nz+1,Nx);      % set the derivative of the nonlinear function to zero 
    Mpsi     = MatVecProd_psi(psi);  
    % compute the residual r(psi) = theta(psi) + M*psi - b
    r(1:Nz,:) = Theta(psi(1:Nz,:)   )    + Mpsi(1:Nz,:)  - rhs(1:Nz,:); 
    r(Nz+1,:) = Height(psi(Nz+1,:),bath) + Mpsi(Nz+1 ,:) - rhs(Nz+1,:); 
    
    r_norm    = sqrt(sum(sum(r.*r))); % outer residual

end