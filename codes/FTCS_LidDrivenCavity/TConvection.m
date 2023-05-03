function T = TConvection(u,v,T,xb)

    global Nx Ny dx dy dt lambda TBC

    fm = zeros(Nx,Ny);
    fp = zeros(Nx,Ny);
    gm = zeros(Nx,Ny);
    gp = zeros(Nx,Ny);

    Tbottom = (TBC*exp(-0.5*(xb-0.5).^2/0.1^2))';

    % x direction, convective part
    fm(1     ,:) = 0; %0.5*u(1   ,:).*( T(1   ,:) + Tlake       ) - 0.5*abs(u(1   ,:)).*( T(1,:)    - Tlake       );
    fm(2:Nx  ,:) = 0.5*u(2:Nx,:).*( T(2:Nx,:) + T(1:Nx-1,:) ) - 0.5*abs(u(2:Nx,:)).*( T(2:Nx,:) - T(1:Nx-1,:) );
    fm(2:Nx  ,:) = fm(2:Nx  ,:) - lambda*(T(2:Nx,:) - T(1:Nx-1,:))/dx;
    fp(1:Nx-1,:) = fm(2:Nx  ,:);
    fp(Nx    ,:) = 0; %0.5*u(Nx+1,:).*( Tlake     + T(Nx    ,:) ) - 0.5*abs(u(Nx+1,:)).*( Tlake     - T(Nx    ,:) );
    
    % y direction
    gm(:,1     ) = 0; %0.5*v(:,1   ).*( T(:,1   ) + Tlake       ) - 0.5*abs(v(:,1   )).*( T(:,1   ) - Tlake       );
    %gm(:,1     ) =-lambda*(T(:,1) - Tbottom)/(dy/2);
    gm(:,2:Ny  ) = 0.5*v(:,2:Ny).*( T(:,2:Ny) + T(:,1:Ny-1) ) - 0.5*abs(v(:,2:Ny)).*( T(:,2:Ny) - T(:,1:Ny-1) );
    gm(:,2:Ny  ) = gm(:,2:Ny  ) - lambda*(T(:,2:Ny) - T(:,1:Ny-1));
    gp(:,1:Ny-1) = gm(:,2:Ny  );
    gp(:,Ny    ) = 0;
    
    % finite volume update
    T = T - dt/dx*( fp - fm ) - dt/dy*( gp - gm );
end
