% Jacobian dQ1/dT of Q1(T)
function jac = dQ1(T)

global hL cS cL rhoS rhoL Tc epsilon

b = T <= Tc - epsilon; % boolean variable
jac = b.* (rhoS*cS);

b = T > Tc - epsilon; % boolean variable
dqdT = (rhoL*cL*((Tc+epsilon) - Tc) + rhoL*hL - rhoS*cS*((Tc-epsilon)-Tc))/(2*epsilon);
jac = jac + b.*dqdT;


end