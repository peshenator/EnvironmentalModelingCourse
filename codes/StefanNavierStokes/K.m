% thermal conductivity as a function of temperature
% in the Stefan problem 
function k = K(T)
global Tc KL KS 

% solid
b = T < Tc;
k = b.*KS;

% liquid
b = T >= Tc;
k = k + b.* KL;
