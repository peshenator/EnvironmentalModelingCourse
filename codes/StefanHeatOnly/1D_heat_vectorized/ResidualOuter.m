function res = ResidualOuter(T,a,b,c,rhs)

global Nx;

res = Q(T) - rhs;
%now add M*T

% res(:) = 0;
res(1 )     = res(1     )                        + b(1 )*T(1 ) + c(1)*T(2);
res(2:Nx-1) = res(2:Nx-1) + a(2:Nx-1).*T(1:Nx-2) + b(2:Nx-1).*T(2:Nx-1) + c(2:Nx-1).*T(3:Nx);
res(Nx)     = res(Nx    ) + a(Nx    ) *T(Nx-1  ) + b(Nx)*T(Nx);
 
end