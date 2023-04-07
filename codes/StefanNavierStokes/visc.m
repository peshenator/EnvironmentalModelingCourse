% a toy model for the temperature dependence of viscosity
function viscosity = visc(T)
global nu Ts epsilonT;
    
    eps = 1e-3;
    viscosity = nu*ones(size(T));
    % nusolid = 100; % solid visocsity
    % b = T >= Ts + eps;
    % viscosity = b.*( nu*exp(abs(b.*T - (Ts-1e-5))^2) );

    % b = T < Ts + eps;
    % viscosity = viscosity + nusolid*exp((Ts-1e-5)^2)*b;
end