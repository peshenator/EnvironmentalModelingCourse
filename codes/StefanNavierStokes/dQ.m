% Jacobian dQ/dT of Q(T)=Q1 - Q2
function jac = dQ(T)

global hL cS cL rhoS rhoL Tc epsilonT

% ice
b = T <= Tc-epsilonT;
jac = b.* (rhoS*cS);
% water
b = T >= Tc+epsilonT;
jac = jac + b.* (rhoL*cL);
% transition regiion
b = (T > Tc - epsilonT) .*  (T < Tc + epsilonT);
dqdT = (rhoL*cL*((Tc+epsilonT) - Tc) + rhoL*hL - rhoS*cS*((Tc-epsilonT)-Tc))/(2*epsilonT);
jac = jac + b.*dqdT;

end