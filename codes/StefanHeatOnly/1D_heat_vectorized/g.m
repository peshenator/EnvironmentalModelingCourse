% function whose root we have to find 
function y=g(gamma)
global kappaL kappaR hs rhoR Ks Kl Ts TL TR

y= gamma*sqrt(kappaL)*hs*rhoR - ...  
   Ks*(Ts-TL)/erf(gamma)*exp(-gamma^2)/sqrt(pi*kappaL) - ... 
   Kl*(Ts-TR)/erfc(gamma*sqrt(kappaL/kappaR))*exp(-gamma^2*kappaL/kappaR)/sqrt(pi*kappaR); 
