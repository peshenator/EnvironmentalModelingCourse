function [Hmx,Hpx,Hmy,Hpy,rhs,Hx,Hy] = LinearPartCoef(Hb,u,v)
global g gamma dt dx dy Nx Ny

Hmx   = zeros(Nx,Ny);  % total water depth on the cell faces [1:Nx]   in x
Hpx   = zeros(Nx,Ny);  % total water depth on the cell faces [2:Nx+1] in x
Hmy   = zeros(Nx,Ny);  % total water depth on the cell faces [1:Ny]   in y
Hpy   = zeros(Nx,Ny);  % total water depth on the cell faces [2:Ny+1] in y

Hmx(1     ,:) = Hb(1,:);
Hmx(2:Nx  ,:) = max(Hb(2:Nx,:) , Hb(1:Nx-1,:)); % total water depth on the cell faces [1:Nx]   in x
Hpx(1:Nx-1,:) = Hmx(2:Nx,:);
Hpx(Nx    ,:) = Hb(Nx,:);

Hmy(:,1     ) = Hb(:,1);
Hmy(:,2:Ny  ) = max(Hb(:,2:Ny) , Hb(:,1:Ny-1)); % total water depth on the cell faces [1:Nx]   in x
Hpy(:,1:Ny-1) = Hmy(:,2:Ny);
Hpy(:,Ny    ) = Hb(:,Ny);

% these we need to store for the velocity update
bool = (Hmx == 0); 
Hx = Hmx./(Hmx + dt*gamma + bool);
bool = (Hmy == 0); 
Hy = Hmy./(Hmy + dt*gamma + bool);

bool = (Hmx == 0);  % this is simply a trick to avoid division by zero when H = 0
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