% Nonlinear volume function 
%   Input: 
%     - eta = free surface elevation 
%     - h   = bottom depth 
%   Output: 
%     - vol = cell volume 
function vol=V(eta,h)
vol = max(0,eta+h);     
