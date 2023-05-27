% auxiliary function that allows to limit velocity in the dry cells
function [u,v] = VelocityFilter(u,v,Hb)

    global Nx Ny;
    
    Hx = zeros(Nx+1,Ny);
    Hy = zeros(Nx,Ny+1);
    
    Hx(1     ,:) = max(Hb(1,:),Hb(Nx,:));   % for periodic BC
    Hx(2:Nx  ,:) = max(Hb(2:Nx,:),Hb(1:Nx-1,:)); % total depth at the cell faces
    Hx(  Nx+1,:) = max(Hb(Nx,:),Hb(1,:));   % for periodic BC
    
    Hy(:,1     ) = Hb(:,1);
    Hy(:,2:Ny  ) = max(Hb(:,2:Ny),Hb(:,1:Ny-1)); % total depth at the cell faces
    Hy(:,  Ny+1) = Hb(:,Ny);

    Hu = Hx.*u;
    Hv = Hy.*v;

    % define iH = 1/H via a filter
    eps = 1e-6;
    iHx = Hx./( Hx.^2 + eps*(Hx + 1) );
    iHy = Hy./( Hy.^2 + eps*(Hy + 1) );

    u = Hu.*iHx;
    v = Hv.*iHy;

end