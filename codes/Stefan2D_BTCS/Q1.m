% define Q1 in Q = Q1-Q2
function q1 = Q1(T)


global hL cS cL rhoS rhoL Tc epsilon;

b = T <= Tc - epsilon; % boolean variable
q1 = b.*( rhoS*cS*(T-Tc) );

% energy in the regularized tranzition zone
% Q(T) = Q(T0) + dQdT(T0)*(T - T0), where T0 = Tc-epsilon
b = T > Tc - epsilon; % boolean variable
dqdT = (rhoL*cL*((Tc+epsilon) - Tc) + rhoL*hL - rhoS*cS*((Tc-epsilon)-Tc))/(2*epsilon);
q1 =q1 + b.*( rhoS*cS*((Tc-epsilon)-Tc) + dqdT*(T - (Tc - epsilon)) );


end