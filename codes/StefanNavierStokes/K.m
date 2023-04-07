% thermal conductivity as a function of temperature
% in the Stefan problem 
function k = K(T)
global Ts KL KS 

% solid
b = T < Ts;
k = b.*KS;

% liquid
b = T >= Ts;
k = k + b.* KL;
