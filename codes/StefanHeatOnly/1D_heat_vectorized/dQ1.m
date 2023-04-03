% Jacobian dQ1/dT of Q1(T)
function jac = dQ1(T)

global hs cL cR rhoL rhoR Ts epsilon

b = T <= Ts - epsilon; % boolean variable
jac = b.* (rhoL*cL);

b = T > Ts - epsilon; % boolean variable
dqdT = (rhoR*cR*((Ts+epsilon) - Ts) + rhoR*hs - rhoL*cL*((Ts-epsilon)-Ts))/(2*epsilon);
jac = jac + b.*dqdT;


end