function y=Thetaf(psi)
global alpha thetas thetar n m Ks psic 

ylogic = (psi<=0);
y = ylogic.*( thetar + (thetas-thetar)./( (1+abs(alpha*psi).^n).^m ) );

ylogic = (psi> 0);
y = y + thetas*ylogic;

end