function y=Theta(psi)
global alpha thetas thetar n m Ks psic 

logic = (psi<=0);
y = logic.*( thetar + (thetas-thetar)./( (1+abs(alpha*psi).^n).^m ) )+...
    thetas*(~logic);

end