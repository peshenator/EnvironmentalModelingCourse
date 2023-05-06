% function whose root we have to find 
function y=g(gamma)
global kappaS kappaL hL rhoL KS KL Tc Tair Tlake

y= gamma*sqrt(kappaS)*hL*rhoL - ...  
   KS*(Tc-Tair)/erf(gamma)*exp(-gamma^2)/sqrt(pi*kappaS) - ... 
   KL*(Tc-Tlake)/erfc(gamma*sqrt(kappaS/kappaL))*exp(-gamma^2*kappaS/kappaL)/sqrt(pi*kappaL); 
