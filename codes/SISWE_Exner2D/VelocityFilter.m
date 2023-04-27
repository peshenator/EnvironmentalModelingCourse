function [u,v] = VelocityFilter(u,v,Hx,Hy)
global Nx Ny;

Hu = Hx.*u;
Hv = Hy.*v;
% define iH = 1/H via a filter in the case where H is close to 0
eps = 1e-6;
iHx = Hx./(Hx.^2 + eps*(Hx + 1));
iHy = Hy./(Hy.^2 + eps*(Hy + 1));

u = Hu.*iHx;
v = Hv.*iHy;
end