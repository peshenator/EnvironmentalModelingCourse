% only convection part of the temperature equation, the diffusion part is implicit 
% and goes to the CGsolver
function T = TConvection(u,v,T,xb)

    global Nx Ny dx dy dt lambda TBC

    fm = zeros(Nx,Ny);
    fp = zeros(Nx,Ny);
    gm = zeros(Nx,Ny);
    gp = zeros(Nx,Ny);

    Tbottom = (TBC*exp(-0.5*(xb-0.5).^2/0.1^2))';

    % x direction
    fm(1     ,:) = 0;   
    fm(2:Nx  ,:) = 0.5*u(2:Nx,:).*( T(2:Nx,:) + T(1:Nx-1,:) ) - 0.5*abs(u(2:Nx,:)).*( T(2:Nx,:) - T(1:Nx-1,:) );
    fp(1:Nx-1,:) = fm(2:Nx  ,:);
    fp(Nx    ,:) = 0;

    % y direction
    gm(:,1     ) = 0;
    % gm(:,1     ) =-lambda*(T(:,1) - Tbottom)/(dy/2);
    gm(:,2:Ny  ) = 0.5*v(:,2:Ny).*( T(:,2:Ny) + T(:,1:Ny-1) ) - 0.5*abs(v(:,2:Ny)).*( T(:,2:Ny) - T(:,1:Ny-1) );
    gp(:,1:Ny-1) = gm(:,2:Ny  );
    gp(:,Ny    ) = 0;
    
    % finite volume update
    T = T - dt/dx*( fp - fm ) - dt/dy*( gp - gm );
end