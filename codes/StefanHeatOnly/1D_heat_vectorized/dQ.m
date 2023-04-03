% Jacobian dQ/dT of Q(T)=Q1 - Q2
function jac = dQ(T)

global hs cL cR rhoL rhoR Ts epsilon

% ice
b = T <= Ts-epsilon;
jac = b.* (rhoL*cL);
% water
b = T >= Ts+epsilon;
jac = jac + b.* (rhoR*cR);
% transition regiion
b = (T > Ts - epsilon) .*  (T < Ts + epsilon);
dqdT = (rhoR*cR*((Ts+epsilon) - Ts) + rhoR*hs - rhoL*cL*((Ts-epsilon)-Ts))/(2*epsilon);
jac = jac + b.*dqdT;

end