% Jacobian dQ1/dT of Q1(T)
function jac = dQ1(T)

global hs cS cL rhoS rhoL Ts epsilonT

b = T <= Ts - epsilonT; % boolean variable
jac = b.* (rhoS*cS);

b = T > Ts - epsilonT; % boolean variable
dqdT = (rhoL*cL*((Ts+epsilonT) - Ts) + rhoL*hs - rhoS*cS*((Ts-epsilonT)-Ts))/(2*epsilonT);
jac = jac + b.*dqdT;


end