% the energy function for the BTCS scheme has to be slightly regularized
% epsilon is a regularizatio parameter
function q = Q(T)

global hL cS cL rhoS rhoL Tc epsilonT;

b = T <= Tc - epsilonT; % boolean variable
q = b.*( rhoS*cS*(T-Tc) );

b = T >= Tc + epsilonT;
q =q + b.*( rhoL*cL*(T-Tc) + rhoL*hL );

b = (T > Tc - epsilonT) .*  (T < Tc + epsilonT);
dqdT = (rhoL*cL*((Tc+epsilonT) - Tc) + rhoL*hL - rhoS*cS*((Tc-epsilonT)-Tc))/(2*epsilonT);
q = q + b.*( rhoS*cS*((Tc-epsilonT)-Tc) + dqdT*(T - (Tc - epsilonT)) );


end