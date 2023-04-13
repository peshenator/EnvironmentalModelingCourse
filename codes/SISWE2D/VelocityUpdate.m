function [u,v] = VelocityUpdate(u,v,etab,Hx,Hy)
    global g dt dx dy Nx Ny;

    kx = g*dt/dx;
    ky = g*dt/dy;

    u(1   ,:) = 0; % left wall
    u(2:Nx,:) = Hx(2:Nx,:).*( u(2:Nx,:) - kx*( etab(2:Nx,:) - etab(1:Nx-1,:) ) );
    u(Nx+1,:) = 0; % right wall

    v(:,1   ) = 0; % left wall
    v(:,2:Ny) = Hy(:,2:Ny).*( v(:,2:Ny) - ky*( etab(:,2:Ny) - etab(:,1:Ny-1) ) );
    v(:,Ny+1) = 0; % right wall

end