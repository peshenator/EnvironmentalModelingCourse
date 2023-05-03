% Jacobian dQ/dT of Q(T)=Q1 - Q2
function jac = dQ(T)

global hs cS cL rhoS rhoL Ts epsilon

% ice
b = T <= Ts-epsilon;
jac = b.* (rhoS*cS);
% water
b = T >= Ts+epsilon;
jac = jac + b.* (rhoL*cL);
% transition regiion
b = (T > Ts - epsilon) .*  (T < Ts + epsilon);
dqdT = (rhoL*cL*((Ts+epsilon) - Ts) + rhoL*hs - rhoS*cS*((Ts-epsilon)-Ts))/(2*epsilon);
jac = jac + b.*dqdT;

end