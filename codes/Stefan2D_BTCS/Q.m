% the energy function for the BTCS scheme has to be slightly regularized
% epsilon is a regularizatio parameter
function q = Q(T)

global hL cS cL rhoS rhoL Tc epsilon;

b = T <= Tc - epsilon; % boolean variable
q = b.*( rhoS*cS*(T-Tc) );

b = T >= Tc + epsilon;
q =q + b.*( rhoL*cL*(T-Tc) + rhoL*hL );

b = (T > Tc - epsilon) .*  (T < Tc + epsilon);
dqdT = (rhoL*cL*((Tc+epsilon) - Tc) + rhoL*hL - rhoS*cS*((Tc-epsilon)-Tc))/(2*epsilon);
q = q + b.*( rhoS*cS*((Tc-epsilon)-Tc) + dqdT*(T - (Tc - epsilon)) );


end