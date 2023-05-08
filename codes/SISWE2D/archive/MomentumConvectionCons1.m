% Explicit subsystem in the conservative formulation
% 
% INPUT : Hub, Hvb, Hb - conservative variables (x- and y-momenta and total water depth)
% OUTPUT:  u,  v - primitive variables (x- and y-velocity)
function [u,v] = MomentumConvectionCons1(Hub,Hvb,Hb,u,v)
global Nx Ny dt dx dy;

F = zeros(Nx+1,Ny,2);
G = zeros(Nx,Ny+1,2);

um = 0.5*(u - abs(u));
up = 0.5*(u + abs(u));
vm = 0.5*(v - abs(v));
vp = 0.5*(v + abs(v));

F(2:Nx,:,1) = um(2:Nx,:).*Hub(2:Nx,:) + up(2:Nx,:).*Hub(1:Nx-1,:);
F(2:Nx,:,2) = um(2:Nx,:).*Hvb(2:Nx,:) + up(2:Nx,:).*Hvb(1:Nx-1,:);

G(:,2:Nx,1) = vm(:,2:Ny).*Hub(:,2:Ny) + vp(:,2:Ny).*Hub(:,1:Ny-1);
G(:,2:Nx,2) = vm(:,2:Ny).*Hvb(:,2:Ny) + vp(:,2:Ny).*Hvb(:,1:Ny-1);
% update the conservative variables in the cell centers
Hub = Hub - dt/dx*( F(2:Nx+1,:,1) - F(1:Nx,:,1) ) - dt/dy*( G(:,2:Ny+1,1) - G(:,1:Ny,1) );
Hvb = Hvb - dt/dx*( F(2:Nx+1,:,2) - F(1:Nx,:,2) ) - dt/dy*( G(:,2:Ny+1,2) - G(:,1:Ny,2) );

% Now, we need to average back the conservative variables to the cell faces
Hu = zeros(Nx+1,Ny);
Hv = zeros(Nx,Ny+1);
Hu(2:Nx,:) = 0.5*( Hub(2:Nx,:) + Hub(1:Nx-1,:) );
Hv(:,2:Ny) = 0.5*( Hvb(:,2:Ny) + Hvb(:,1:Ny-1) );

Hx = zeros(Nx+1,Ny);
Hy = zeros(Nx,Ny+1);
Hx(2:Nx,:) = 0.5*( Hb(2:Nx,:) + Hb(1:Nx-1,:) );
Hx(1   ,:) = Hx(2 ,:);
Hx(Nx+1,:) = Hx(Nx,:);

Hy(:,2:Ny) = 0.5*( Hb(:,2:Ny) + Hb(:,1:Ny-1) );
Hy(:,1   ) = Hy(:,2 );
Hy(:,Ny+1) = Hy(:,Ny);
% define iH = 1/H via a filter in the case where H is close to 0
eps = 1e-6;
iHx = Hx./(Hx.^2 + eps*(Hx + 1));
iHy = Hy./(Hy.^2 + eps*(Hy + 1));

u = Hu.*iHx;
v = Hv.*iHy;
end