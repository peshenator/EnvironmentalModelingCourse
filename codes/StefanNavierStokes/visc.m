% a toy model for the temperature dependence of viscosity
function viscosity = visc(T)
global nu Tc epsilonT;
    
%     eps = 1e-3;
%     viscosity = nu*ones(size(T));
    nusolid = 1000; % solid visocsity
    b = T >= Tc + epsilonT;
    viscosity = b.*( nu*exp(abs(b.*T - (Tc-1e-5))^2) );

    b = T < Tc + epsilonT;
    viscosity = viscosity + nusolid*exp((Tc-1e-5)^2)*b;
end