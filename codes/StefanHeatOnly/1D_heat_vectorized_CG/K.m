% thermal conductivity as a function of temperature
% in the Stefan problem 
function k = K(T)
global Ts Kl Ks 

% solid
b = T < Ts;
k = b.*Ks;

% liquid
b = T >= Ts;
k = k + b.* Kl;
