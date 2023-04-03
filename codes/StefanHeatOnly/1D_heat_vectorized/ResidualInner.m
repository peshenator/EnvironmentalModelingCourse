function [res,dQ12] = ResidualInner(T,Talpha,a,b,c,rhs)

global Nx;

res = Q1(T) - Q2(Talpha) - dQ2(Talpha).*(T - Talpha) - rhs;
% add M*T
res(1) = res(1) + b(1)*T(1) + c(1)*T(2);
res(Nx) = res(Nx) + a(Nx)*T(Nx-1) + b(Nx)*T(Nx);
res(2:Nx-1) = res(2:Nx-1) + a(2:Nx-1).*T(1:Nx-2) + b(2:Nx-1).*T(2:Nx-1) + c(2:Nx-1).*T(3:Nx);

dQ12 = dQ1(T) - dQ2(Talpha);

end