% define Q1 in Q = Q1-Q2
function q1 = Q1(T)


global hL cS cL rhoS rhoL Ts epsilon;

b = T <= Ts - epsilon; % boolean variable
q1 = b.*( rhoS*cS*(T-Ts) );

% energy in the regularized tranzition zone
% Q(T) = Q(T0) + dQdT(T0)*(T - T0), where T0 = Ts-epsilon
b = T > Ts - epsilon; % boolean variable
dqdT = (rhoL*cL*((Ts+epsilon) - Ts) + rhoL*hL - rhoS*cS*((Ts-epsilon)-Ts))/(2*epsilon);
q1 =q1 + b.*( rhoS*cS*((Ts-epsilon)-Ts) + dqdT*(T - (Ts - epsilon)) );


end