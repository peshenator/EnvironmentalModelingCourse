% define Q1 in Q = Q1-Q2
function q1 = Q1(T)


global hs cS cL rhoS rhoL Ts epsilonT;

b = T <= Ts - epsilonT; % boolean variable
q1 = b.*( rhoS*cS*(T-Ts) );

% energy in the regularized tranzition zone
% Q(T) = Q(T0) + dQdT(T0)*(T - T0), where T0 = Ts-epsilon
b = T > Ts - epsilonT; % boolean variable
dqdT = (rhoL*cL*((Ts+epsilonT) - Ts) + rhoL*hs - rhoS*cS*((Ts-epsilonT)-Ts))/(2*epsilonT);
q1 =q1 + b.*( rhoS*cS*((Ts-epsilonT)-Ts) + dqdT*(T - (Ts - epsilonT)) );


end