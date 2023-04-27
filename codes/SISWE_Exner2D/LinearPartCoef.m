function [Hmx,Hpx,Hmy,Hpy,Hx,Hy,tHx,tHy,rhs] = LinearPartCoef(Hb,u,v)
global g gamma dt dx dy Nx Ny

Hmx   = zeros(Nx,Ny);  % total water depth on the cell faces [1:Nx]   in x
Hpx   = zeros(Nx,Ny);  % total water depth on the cell faces [2:Nx+1] in x
Hmy   = zeros(Nx,Ny);  % total water depth on the cell faces [1:Ny]   in y
Hpy   = zeros(Nx,Ny);  % total water depth on the cell faces [2:Ny+1] in y

Hmx(1     ,:) = max(Hb(1   ,:) , Hb(  Nx  ,:));
Hmx(2:Nx  ,:) = max(Hb(2:Nx,:) , Hb(1:Nx-1,:)); % total water depth on the cell faces [1:Nx]   in x
Hpx(1:Nx-1,:) = Hmx(2:Nx,:);
Hpx(Nx    ,:) = max(Hb(Nx  ,:) , Hb(1     ,:));

tHx = zeros(Nx+1,Ny);
tHx(1:Nx,:) = Hmx;
tHx(Nx+1,:) = Hpx(Nx,:);

Hmy(:,1     ) = Hb(:,1);
Hmy(:,2:Ny  ) = max(Hb(:,2:Ny) , Hb(:,1:Ny-1)); % total water depth on the cell faces [1:Nx]   in x
Hpy(:,1:Ny-1) = Hmy(:,2:Ny);
Hpy(:,Ny    ) = Hb(:,Ny);

tHy = zeros(Nx,Ny+1);
tHy(:,1:Ny) = Hmy;
tHy(:,Ny+1) = Hpy(:,Ny);

% these we need to store for the velocity update
Hx = zeros(Nx+1,Ny);
Hy = zeros(Nx,Ny+1);
bool = (Hmx == 0); 
Hx(1:Nx,:) = Hmx./(Hmx + dt*gamma + bool);
Hx(Nx+1,:) = Hpx(Nx,:)./(Hpx(Nx,:) + dt*gamma + (Hpx(Nx,:) == 0));
bool = (Hmy == 0); 
Hy(:,1:Ny) = Hmy./(Hmy + dt*gamma + bool);
Hy(:,Ny+1) = Hpy(:,Ny)./(Hpy(:,Ny) + dt*gamma + (Hpy(:,Ny) == 0));

bool = (Hmx == 0); 
Hmx = Hmx.^2./(Hmx + dt*gamma + bool);
bool = (Hpx == 0); 
Hpx = Hpx.^2./(Hpx + dt*gamma + bool);
bool = (Hmy == 0); 
Hmy = Hmy.^2./(Hmy + dt*gamma + bool);
bool = (Hpy == 0); 
Hpy = Hpy.^2./(Hpy + dt*gamma + bool);

% compute the right hand side of the eta system

% inner points
dtdx = dt/dx;
dtdy = dt/dy;
rhs = Hb  - dtdx*( Hpx.*u(2:Nx+1,:) - Hmx.*u(1:Nx,:) ) - ...
            dtdy*( Hpy.*v(:,2:Ny+1) - Hmy.*v(:,1:Ny) );
end




% Hmx(1     ,:) = max(0,etab(1,:) + bathx(1,:));
% Hmx(2:Nx  ,:) = max(0,0.5*( etab(2:Nx,:) + etab(1:Nx-1,:) ) + bathx(2:Nx,:)); % total water depth on the cell faces [1:Nx]   in x
% Hpx(1:Nx-1,:) = Hmx(2:Nx,:);
% Hpx(Nx    ,:) = max(0,etab(Nx,:) + bathx(Nx+1,:));
% 
% Hmy(:,1     ) = max(0,etab(:,1) + bathy(:,1));
% Hmy(:,2:Ny  ) = max(0,0.5*( etab(:,2:Ny) + etab(:,1:Ny-1) ) + bathy(:,2:Ny)); % total water depth on the cell faces [1:Nx]   in x
% Hpy(:,1:Ny-1) = Hmx(:,2:Ny);
% Hpy(:,Ny    ) = max(0,etab(:,Ny) + bathy(:,Ny+1));