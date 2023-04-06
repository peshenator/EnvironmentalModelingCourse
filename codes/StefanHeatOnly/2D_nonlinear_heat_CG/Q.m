% the energy function for the BTCS scheme has to be slightly regularized
% epsilon is a regularizatio parameter
function q = Q(T)

global hs cS cL rhoS rhoL Ts epsilon;

b = T <= Ts - epsilon; % boolean variable
q = b.*( rhoS*cS*(T-Ts) );

b = T >= Ts + epsilon;
q =q + b.*( rhoL*cL*(T-Ts) + rhoL*hs );

b = (T > Ts - epsilon) .*  (T < Ts + epsilon);
dqdT = (rhoL*cL*((Ts+epsilon) - Ts) + rhoL*hs - rhoS*cS*((Ts-epsilon)-Ts))/(2*epsilon);
q = q + b.*( rhoS*cS*((Ts-epsilon)-Ts) + dqdT*(T - (Ts - epsilon)) );


end