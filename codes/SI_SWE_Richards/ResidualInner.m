% compute the residual of the inner iterations and dC12 = dC1^{beta+1} - dC2^alpha
function [r,r_norm] = ResidualInner(psi,psik,bath,rhs)
global Nx Nz Kx Kz di

    r = zeros(Nz+1,Nx);     % residual of the linearized system
    di   = zeros(Nz+1,Nx);  
    Mpsi = MatVecProd_psi(psi);   % linear part of the system 
    
    r(1:Nz,:) = Theta1(psi(1:Nz,:)) - (Theta2(psik(1:Nz,:)) + dTheta2(psik(1:Nz,:)).*(psi(1:Nz,:)-psik(1:Nz,:))) + Mpsi(1:Nz,:) - rhs(1:Nz,:);  
    di(1:Nz,:) = dTheta1(psi(1:Nz,:)) - dTheta2(psik(1:Nz,:)); 
    
    r(Nz+1,:) = Height(psi(Nz+1,:),bath) + Mpsi(Nz+1,:) - rhs(Nz+1,:); % -(V2(psik(k,i))+dV2(psik(k,i))*(psi(k,i)-psik(k,i)))
    di(Nz+1,:) = dHeight(psi(Nz+1,:),bath);                      % -dTheta2(psik(k,i))

    for i=1:Nx
        for k=1:Nz+1
            if(Kx(k,i) + Kx(k,i+1) + Kz(k,i) + Kz(k+1,i) == 0)
                r(k ,i) = 0;
                di(k,i) = 1; 
            end
        end
    end

    r_norm = sqrt(sum(sum(r.*r))); 

end