% Jacobian dQ/dT of Q(T)=Q1 - Q2
function jac = dQ(T)

global hL cS cL rhoS rhoL Tc epsilon

% ice
b = T <= Tc-epsilon;
jac = b.* (rhoS*cS);
% water
b = T >= Tc+epsilon;
jac = jac + b.* (rhoL*cL);
% transition regiion
b = (T > Tc - epsilon) .*  (T < Tc + epsilon);
dqdT = (rhoL*cL*((Tc+epsilon) - Tc) + rhoL*hL - rhoS*cS*((Tc-epsilon)-Tc))/(2*epsilon);
jac = jac + b.*dqdT;

end