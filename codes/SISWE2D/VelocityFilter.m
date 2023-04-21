function [u,v] = VelocityFilter(u,v,Hb)
global Nx Ny;

% Now, we need to average back the conservative variables to the cell faces

Hx = zeros(Nx+1,Ny);
Hy = zeros(Nx,Ny+1);
Hx(2:Nx,:) = 0.5*( Hb(2:Nx,:) + Hb(1:Nx-1,:) );
Hx(1   ,:) = Hx(2 ,:);
Hx(Nx+1,:) = Hx(Nx,:);

Hy(:,2:Ny) = 0.5*( Hb(:,2:Ny) + Hb(:,1:Ny-1) );
Hy(:,1   ) = Hy(:,2 );
Hy(:,Ny+1) = Hy(:,Ny);

Hu = Hx.*u;
Hv = Hy.*v;
% define iH = 1/H via a filter in the case where H is close to 0
eps = 1e-6;
iHx = Hx./(Hx.^2 + eps*(Hx + 1));
iHy = Hy./(Hy.^2 + eps*(Hy + 1));

u = Hu.*iHx;
v = Hv.*iHy;
end