function T = TConvection(u,v,T)

    global Nx Ny dx dy dt Tlake Tair

    fm = zeros(Nx,Ny);
    fp = zeros(Nx,Ny);
    gm = zeros(Nx,Ny);
    gp = zeros(Nx,Ny);

    % x direction
    fm(1     ,:) = 0; %0.5*u(1   ,:).*( T(1   ,:) + Tlake       ) - 0.5*abs(u(1   ,:)).*( T(1,:)    - Tlake       );
    fm(2:Nx  ,:) = 0.5*u(2:Nx,:).*( T(2:Nx,:) + T(1:Nx-1,:) ) - 0.5*abs(u(2:Nx,:)).*( T(2:Nx,:) - T(1:Nx-1,:) );
    fp(1:Nx-1,:) = 0.5*u(2:Nx,:).*( T(2:Nx,:) + T(1:Nx-1,:) ) - 0.5*abs(u(2:Nx,:)).*( T(2:Nx,:) - T(1:Nx-1,:) );
    fp(Nx    ,:) = 0; %0.5*u(Nx+1,:).*( Tlake     + T(Nx    ,:) ) - 0.5*abs(u(Nx+1,:)).*( Tlake     - T(Nx    ,:) );

    % y direction
    gm(:,1     ) = 0; %0.5*v(:,1   ).*( T(:,1   ) + Tlake       ) - 0.5*abs(v(:,1   )).*( T(:,1   ) - Tlake       );
    gm(:,2:Ny  ) = 0.5*v(:,2:Ny).*( T(:,2:Ny) + T(:,1:Ny-1) ) - 0.5*abs(v(:,2:Ny)).*( T(:,2:Ny) - T(:,1:Ny-1) );
    gp(:,1:Ny-1) = 0.5*v(:,2:Ny).*( T(:,2:Ny) + T(:,1:Ny-1) ) - 0.5*abs(v(:,2:Ny)).*( T(:,2:Ny) - T(:,1:Ny-1) );
    gp(:,Ny    ) = 0.5*v(:,Ny+1).*( Tair      + T(:,Ny    ) ) - 0.5*abs(v(:,Ny+1)).*( Tair      - T(:,Ny    ) );
    
    % finite volume update
    T = T - dt/dx*( fp - fm ) - dt/dy*( gp - gm );
end
