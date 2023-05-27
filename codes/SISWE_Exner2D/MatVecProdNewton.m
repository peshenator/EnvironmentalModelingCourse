function Meta = MatVecProdNewton(etab,bathb,Hxm,Hxp,Hym,Hyp,wet)

    global g dt dx dy Nx Ny;

    kx = g*(dt/dx)^2;
    ky = g*(dt/dy)^2;
    
    fm = zeros(Nx,Ny);
    fp = zeros(Nx,Ny);
    
    % fluxes in the X directin
    fm(1     ,:) = Hxm(1     ,:).*( etab(1   ,:) - etab(Nx    ,:) );    % periodic BC
    fm(2:Nx  ,:) = Hxm(2:Nx  ,:).*( etab(2:Nx,:) - etab(1:Nx-1,:) );
    fp(1:Nx-1,:) = fm(2:Nx,:);
    fp(Nx    ,:) = Hxp( Nx   ,:).*( etab(1   ,:) - etab(Nx    ,:) );    % periodic BC

    Meta =-kx*( fp - fm );

    % fluxes in the Y directin
    fm(:,1     ) = 0;
    fm(:,2:Ny  ) = Hym(:,2:Ny).*( etab(:,2:Ny) - etab(:,1:Ny-1) );
    fp(:,1:Ny-1) = fm(:,2:Ny);
    fp(:,  Ny  ) = 0;

    Meta = Meta - ky*( fp - fm );

    % for the Newton iterations we need to add p(eta)*eta to of M*eta, see
    % the Matrix of the Newton method in the lecture
    Meta = Meta + wet.*etab;
    
end