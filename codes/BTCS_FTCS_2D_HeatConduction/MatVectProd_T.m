% matrix-vector product for the 2D linear heat equation
function AT = MatVectProd_T(T,Lxm,Lxp,Lym,Lyp)
    global Nx Ny dt dx dy;
    global BCDirichletXL BCDirichletXR BCDirichletYL BCDirichletYR;

    fm = zeros(Nx,Ny);
    fp = zeros(Nx,Ny);
    gm = zeros(Nx,Ny);
    gp = zeros(Nx,Ny);

    % numerical flux in the x-direction
    fm(2:Nx  ,:) =-Lxm(2:Nx,:).*( T(2:Nx,:) - T(1:Nx-1,:) )/dx;   % Fourier law
    fp(1:Nx-1,:) = fm(2:Nx,:);
    
    % numerical flux in the y-direction
    gm(:,2:Ny  ) =-Lym(:,2:Nx).*( T(:,2:Ny) - T(:,1:Ny-1) )/dy;   % Fourier law
    gp(:,1:Ny-1) = gm(:,2:Ny);

    % If Neumann type BC then just set the fluxes fm(1,:), fp(Nx,:), etc to a certain value, we keep it zero here

    % Dirichlet Boundary conditions :
    if (BCDirichletXL) 
        fm(1     ,:) =-Lxm(1 ,:).*T(1 ,:)/dx;  % boundary condition, no heat flux at the left boundary
    end
    if (BCDirichletXR)
        fp(Nx    ,:) = Lxp(Nx,:).*T(Nx,:)/dx;  % boundary condition, no heat flux at the right boundary
    end
    if (BCDirichletYL)
        gm(:,1     ) =-Lym(:,1 ).*T(:,1 )/dy;  % boundary condition, no heat flux at the left boundary
    end
    if (BCDirichletYR)
        gp(:,  Ny  ) = Lym(:,Ny).*T(:,Ny)/dy;  % boundary condition, no heat flux at the right boundary
    end
    
    AT = T  + dt/dx*( fp - fm ) + dt/dy*( gp - gm );
end