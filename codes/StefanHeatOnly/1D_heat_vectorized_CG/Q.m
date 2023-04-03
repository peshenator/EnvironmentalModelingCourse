% the energy function for the BTCS scheme has to be slightly regulzrized
function q = Q(T)

global hs cL cR rhoL rhoR Ts epsilon;

b = T <= Ts - epsilon; % boolean variable
q = b.*( rhoL*cL*(T-Ts) );

b = T >= Ts + epsilon;
q =q + b.*( rhoR*cR*(T-Ts) + rhoR*hs );

b = (T > Ts - epsilon) .*  (T < Ts + epsilon);
dqdT = (rhoR*cR*((Ts+epsilon) - Ts) + rhoR*hs - rhoL*cL*((Ts-epsilon)-Ts))/(2*epsilon);
q = q + b.*( rhoL*cL*((Ts-epsilon)-Ts) + dqdT*(T - (Ts - epsilon)) );


end