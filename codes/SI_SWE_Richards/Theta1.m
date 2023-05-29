function y = Theta1(psi)
global alpha thetas thetar n m Ks psic 


logic = (psi <= psic);
y = Theta(psi).*logic + (Theta(psic) + dTheta(psic)*(psi-psic)).*(~logic);

end

% if(psi<=psic)
%     y = Thetaf(psi); 
% else
%     y = Thetaf(psic) + dTheta(psic)*(psi-psic); 
% end