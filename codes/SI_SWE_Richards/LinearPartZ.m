function [a,b,c,rhs] = LinearPartZ(theta,Fu,Hx,Kz)
global dt dx dz Nx Nz

a   = zeros(Nz+1,Nx);
b   = zeros(Nz+1,Nx);
c   = zeros(Nz+1,Nx);
rhs = zeros(Nz+1,Nx);
dtdx = dt/dx;
dtdz = dt/dz;
dtdz2= dt/dz^2;

% For the free surface equation (k==Nz+1), there is a zero mass flux on the upper boundary = Neumann BC, but there 
% is mass flux through the bottom 

a(1   ,:) = 0;
a(2:Nz,:) =-dtdz2*Kz(2:Nz,:);
a(Nz+1,:) =-dtdz2*Kz(Nz+1,:);

b(1   ,:) =+dtdz2*(Kz(2,:) + 0);
b(2:Nz,:) =+dtdz2*(Kz(3:Nz+1,:) + Kz(2:Nz,:));
b(Nz+1,:) =+dtdz2*(0 + Kz(Nz+1,:));

c(1   ,:) =-dtdz2*Kz(2,:);
c(2:Nz,:) =-dtdz2*Kz(3:Nz+1,:);
c(Nz+1,:) = 0;

rhs(1   ,:) = theta(1   ,:) + dtdz*(Kz(2,:)      - 0         );
rhs(2:Nz,:) = theta(2:Nz,:) + dtdz*(Kz(3:Nz+1,:) - Kz(2:Nz,:)); 
rhs(Nz+1,:) = theta(Nz+1,:) + dtdz*(Kz(Nz+2  ,:) - Kz(Nz+1,:)); 
% x fluxes with homogeneous Neumann BC 
% = do nothing for Richards, but do something for the shallow water part ! 
rhs(Nz+1,:) = rhs(Nz+1,:) - dtdx*(Hx(2:Nx+1).*Fu(2:Nx+1)-Hx(1:Nx).*Fu(1:Nx)); 

end