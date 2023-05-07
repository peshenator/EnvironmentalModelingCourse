% define Q1 in Q = Q1-Q2
function q1 = Q1(T)


global hL cS cL rhoS rhoL Tc epsilonT;

b = T <= Tc - epsilonT; % boolean variable
q1 = b.*( rhoS*cS*(T-Tc) );

% energy in the regularized tranzition zone
% Q(T) = Q(T0) + dQdT(T0)*(T - T0), where T0 = Tc-epsilon
b = T > Tc - epsilonT; % boolean variable
dqdT = (rhoL*cL*((Tc+epsilonT) - Tc) + rhoL*hL - rhoS*cS*((Tc-epsilonT)-Tc))/(2*epsilonT);
q1 =q1 + b.*( rhoS*cS*((Tc-epsilonT)-Tc) + dqdT*(T - (Tc - epsilonT)) );


end