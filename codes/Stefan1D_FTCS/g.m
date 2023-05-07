% function whose root we have to find 
function y=g(gamma)
global lambdaS lambdaL hL rhoL Tc TS TL KS KL;

y= gamma*sqrt(lambdaS)*hL*rhoL - ...  
   KS*(Tc-TS)/erf(gamma)*exp(-gamma^2)/sqrt(pi*lambdaS) - ... 
   KL*(Tc-TL)/erfc(gamma*sqrt(lambdaS/lambdaL))*exp(-gamma^2*lambdaS/lambdaL)/sqrt(pi*lambdaL); 
