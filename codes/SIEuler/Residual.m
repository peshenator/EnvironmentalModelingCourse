function res = Residual(alpha,beta,gamma,d,sol)
global Nx;

residual = zeros(1,Nx);

residual(2:Nx-1) = alpha(2:Nx-1).*sol(1:Nx-2) + beta(2:Nx-1).*sol(2:Nx-1) ...
                   + gamma(2:Nx-1).*sol(3:Nx) - d(2:Nx-1);
residual(1) = beta(1)*sol(1) + gamma(1)*sol(1+1) - d(1);
residual(Nx) = alpha(Nx)*sol(Nx-1) + beta(Nx)*sol(Nx)      - d(Nx);


norm_sol = norm(sol);
if (norm_sol == 0)
    res = 0;
else
    res = norm(residual)/norm_sol;
end

end