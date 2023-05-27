% update the velocity after the Neton's iteration on eta^{n+1}
% here, Hx is H with "^" in the lecture notations on x-faces
% and Hy is H with "^" in the lecture notations on y-faces
function  [u,v] = VelocityUpdate(u,v,etab,Hx,Hy)

    global g dt dx dy Nx Ny;

    u(1   ,:) = Hx(1   ,:).*( u(1   ,:) - g*dt/dx*( etab(1   ,:) - etab(  Nx  ,:) ) );  % periodic BC
    u(2:Nx,:) = Hx(2:Nx,:).*( u(2:Nx,:) - g*dt/dx*( etab(2:Nx,:) - etab(1:Nx-1,:) ) );
    u(Nx+1,:) = Hx(Nx+1,:).*( u(Nx+1,:) - g*dt/dx*( etab(1   ,:) - etab(  Nx  ,:) ) );  % periodic BC

    v(:,1   ) = 0;  % wall BC
    v(:,2:Ny) = Hy(:,2:Ny).*( v(:,2:Ny) - g*dt/dy*( etab(:,2:Ny) - etab(:,1:Ny-1) ) );
    v(:,Ny+1) = 0;  % wall BC

end