% Jacobian dQ1/dT of Q1(T)
function jac = dQ1(T)

global hL cS cL rhoS rhoL Tc epsilonT

b = T <= Tc - epsilonT; % boolean variable
jac = b.* (rhoS*cS);

b = T > Tc - epsilonT; % boolean variable
dqdT = (rhoL*cL*((Tc+epsilonT) - Tc) + rhoL*hL - rhoS*cS*((Tc-epsilonT)-Tc))/(2*epsilonT);
jac = jac + b.*dqdT;


end