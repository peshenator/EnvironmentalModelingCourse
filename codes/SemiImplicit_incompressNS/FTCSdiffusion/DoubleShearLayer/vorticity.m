function omegaz = vorticity(u,v,Nx,Ny,dx,dy)

omegaz = zeros(Nx+1,Ny+1);
% iner vortex points of the domain
omegaz(2:Nx,2:Ny) = (v(2:Nx,2:Ny) - v(1:Nx-1,2:Ny))/dx - (u(2:Nx,2:Ny) - u(2:Nx,1:Ny-1))/dy;

% x left
omegaz(1,2:Ny)    = (v(1,2:Ny) - v(Nx,2:Ny))/dx - (u(1,2:Ny) - u(1,1:Ny-1))/dy;
% x right
omegaz(Nx+1,2:Ny) = (v(1,2:Ny) - v(Nx,2:Ny))/dx - (u(Nx+1,2:Ny) - u(Nx+1,1:Ny-1))/dy;

% ybottom
omegaz(2:Nx,1) = (v(2:Nx,1) - v(1:Nx-1,1))/dx - (u(2:Nx,1) - u(2:Nx,Ny))/dy;
% y top
omegaz(2:Nx,Ny+1) = (v(2:Nx,Ny+1) - v(1:Nx-1,Ny+1))/dx - (u(2:Nx,1) - u(2:Nx,Ny))/dy;

% corners of the domain
omegaz(1,1)    = (v(1,1) - v(Nx,1))/dx - (u(1,1) - u(Nx,1))/dy;
omegaz(Nx+1,1) = (v(1,1) - v(Nx,1))/dx - (u(Nx+1,1) - u(Nx+1,Ny))/dy;

omegaz(Nx+1,Ny+1) = (v(1,Ny+1) - v(Nx,Ny+1))/dx - (u(Nx+1,1) - u(Nx+1,Ny))/dy;
omegaz(1   ,Ny+1) = (v(1,Ny+1) - v(Nx,Ny+1))/dx - (u(1,1) - u(Nx,1))/dy;

end