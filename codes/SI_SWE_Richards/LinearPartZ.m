function [a,b,c,theta,rhs] = LinearPartZ(psi,theta,bath,Fu)
global dt dx dz Nx Nz Kx Kz Hx K g gamma

    a   = zeros(Nz+1,Nx);
    b   = zeros(Nz+1,Nx);
    c   = zeros(Nz+1,Nx);
    rhs = zeros(Nz+1,Nx);
    dtdx = dt/dx;
    dtdz = dt/dz;
    dtdz2= dt/dz^2;

    % Compute volume and hydraulic conductivity at the cell-centers
    theta(1:Nz,:) = Theta(psi(1:Nz,:));         % porous medium
    theta(Nz+1,:) = Height(psi(Nz+1,:),bath);   % free surface layer
    K(1:Nz,:) = Kfun(psi(1:Nz,:));              % porous medium
    K(Nz+1,:) = Kfun(psi(Nz+1,:));              % free surface layer

    % Compute depth of the free-surface layer at the cell-interfaces:
    H = theta(Nz+1,:);
    Hx(1   ) = 0;
    Hx(2:Nx) = max(0,max(H(2:Nx),H(1:Nx-1)));
    Hx(Nx+1) = 0;

    % Compute hydraulic conductivities at the cell-interfaces 
    % inside of the porous medium:
    Kx(1:Nz,1   ) = 0;
    Kx(1:Nz,2:Nx) = max(K(1:Nz,2:Nx),K(1:Nz,1:Nx-1));
    Kx(1:Nz,Nx+1) = 0;
    % artificial hydraulic conductivity in the free surface:
    Kx(Nz+1,:)    = g*dt*Hx.^2./(Hx + dt*gamma + 1e-14);

    Kz(1     ,1:Nx) = 0;
    Kz(2:Nz+1,1:Nx) = max(K(2:Nz+1,1:Nx),K(1:Nz,1:Nx));
    Kz(Nz+2  ,1:Nx) = 0;
    %Kz(Nz+1,:)=0e-13;

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