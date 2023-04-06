% a toy model for the temperature dependence of viscosity
function viscosity = visc(T)
global nu;
    
    viscosity = nu*ones(size(T));    %*exp(1./T);

end