function [M,RHS] = PoissonMatrixPeriodicBC(P,rhs,Nx,Ny,dx2,dy2)

    Nxy = Nx*Ny;
    M = zeros(Nxy,Nxy);
    dP = P;
    eps = 1e-3;
    for k = 1:Nxy
        dP(k) = P(k) + eps;
        M(:,k) = (PoissonEqPressure(dP,Nx,Ny,dx2,dy2) - PoissonEqPressure(P,Nx,Ny,dx2,dy2))/eps;
        dP(k) = dP(k) - eps;
    end
    RHS = reshape(rhs,Nxy,1);
end