% Jacobian dQ1/dT of Q1(T)
function jac = dQ1(T)

global hs cS cL rhoS rhoL Ts epsilon

b = T <= Ts - epsilon; % boolean variable
jac = b.* (rhoS*cS);

b = T > Ts - epsilon; % boolean variable
dqdT = (rhoL*cL*((Ts+epsilon) - Ts) + rhoL*hs - rhoS*cS*((Ts-epsilon)-Ts))/(2*epsilon);
jac = jac + b.*dqdT;


end