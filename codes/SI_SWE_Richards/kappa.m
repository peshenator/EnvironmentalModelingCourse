function K=kappa(psi)
global alpha thetas thetar n m Ks psic 
sat = (Thetaf(psi)-thetar)/(thetas-thetar); % saturation 
K = Ks*sqrt(sat).*(1-(1-sat.^(1/m)).^m).^2;
