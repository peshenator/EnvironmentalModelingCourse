% derivative of the nonlinear volume function 
% Input: free surface elevation (eta) and bottom depth (h)
% Output: derivative of the volume (dvol) 
function dvol = dV(eta,h) 

dvol = 1.0*(eta+h>0);

end