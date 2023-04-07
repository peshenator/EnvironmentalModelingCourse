% the energy function for the BTCS scheme has to be slightly regularized
% epsilon is a regularizatio parameter
function q = Q(T)

global hs cS cL rhoS rhoL Ts epsilonT;

b = T <= Ts - epsilonT; % boolean variable
q = b.*( rhoS*cS*(T-Ts) );

b = T >= Ts + epsilonT;
q =q + b.*( rhoL*cL*(T-Ts) + rhoL*hs );

b = (T > Ts - epsilonT) .*  (T < Ts + epsilonT);
dqdT = (rhoL*cL*((Ts+epsilonT) - Ts) + rhoL*hs - rhoS*cS*((Ts-epsilonT)-Ts))/(2*epsilonT);
q = q + b.*( rhoS*cS*((Ts-epsilonT)-Ts) + dqdT*(T - (Ts - epsilonT)) );


end