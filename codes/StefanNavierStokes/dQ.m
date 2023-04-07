% Jacobian dQ/dT of Q(T)=Q1 - Q2
function jac = dQ(T)

global hs cS cL rhoS rhoL Ts epsilonT

% ice
b = T <= Ts-epsilonT;
jac = b.* (rhoS*cS);
% water
b = T >= Ts+epsilonT;
jac = jac + b.* (rhoL*cL);
% transition regiion
b = (T > Ts - epsilonT) .*  (T < Ts + epsilonT);
dqdT = (rhoL*cL*((Ts+epsilonT) - Ts) + rhoL*hs - rhoS*cS*((Ts-epsilonT)-Ts))/(2*epsilonT);
jac = jac + b.*dqdT;

end