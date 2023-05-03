% function whose root we have to find 
function y=g(gamma)
global kappaS kappaL hs rhoL Tc TS TL KS KL;

y= gamma*sqrt(kappaS)*hs*rhoL - ...  
   KS*(Tc-TS)/erf(gamma)*exp(-gamma^2)/sqrt(pi*kappaS) - ... 
   KL*(Tc-TL)/erfc(gamma*sqrt(kappaS/kappaL))*exp(-gamma^2*kappaS/kappaL)/sqrt(pi*kappaL); 
