function [theta,rhs] = LinearPartCoeff(psi,theta,bath,ustar)
global dt dx dz Nx Nz Kx Kz Hx K g gamma

    rhs = zeros(Nz+1,Nx);
    dtdx = dt/dx;
    dtdz = dt/dz;

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
    Hx_tilde = Hx.^2./(Hx + dt*gamma + 1e-14);  % H_tilde in the lecture notations

    % Compute hydraulic conductivities at the cell-interfaces 
    % inside of the porous medium:
    Kx(1:Nz,1   ) = 0;
    Kx(1:Nz,2:Nx) = max(K(1:Nz,2:Nx),K(1:Nz,1:Nx-1));
    Kx(1:Nz,Nx+1) = 0;
    % artificial hydraulic conductivity in the free surface:
    Kx(Nz+1,:)    = g*dt*Hx_tilde;

    Kz(1     ,1:Nx) = 0;
    Kz(2:Nz+1,1:Nx) = max(K(2:Nz+1,1:Nx),K(1:Nz,1:Nx));
    Kz(Nz+2  ,1:Nx) = 0;
    %Kz(Nz+1,:)=0e-13;
   
    rhs(1   ,:) = theta(1   ,:) + dtdz*(Kz(2,:)      - 0         );
    rhs(2:Nz,:) = theta(2:Nz,:) + dtdz*(Kz(3:Nz+1,:) - Kz(2:Nz,:)); 
    rhs(Nz+1,:) = theta(Nz+1,:) + dtdz*(Kz(Nz+2  ,:) - Kz(Nz+1,:)); 
    % x fluxes with homogeneous Neumann BC 
    % = do nothing for Richards, but do something for the shallow water part ! 
    rhs(Nz+1,:) = rhs(Nz+1,:) - dtdx*(Hx_tilde(2:Nx+1).*ustar(2:Nx+1) - Hx_tilde(1:Nx).*ustar(1:Nx)); 
    
end