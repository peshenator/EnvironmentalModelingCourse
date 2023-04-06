% Result of the matrix-vector product for the Poisson equation

function [px,py]=gradP(p)
global Nx Ny dx dy dt

px = zeros(Nx+1,Ny  );
py = zeros(Ny  ,Ny+1);

% x-direction 
px(1,:)      = 0; %(p(1,:)     - 0 )/(dx/2);    % Left BC
px(Nx+1,:) = 0; %(p(Nx,:)  - 0 )/(dx/2);    % Right BC
px(2:Nx,:) = (dt/dx)*(p(2:Nx,:)- p(1:Nx-1,:));

% y-direction 
py(:,1)      = 0; %(p(:,1)     - 0 )/(dy/2);    % Bottom BC
py(:,Ny+1) = 0; %(p(:,Ny)  - 0 )/(dy/2);    % Top BC
py(:,2:Ny) = (dt/dy)*(p(:,2:Ny)- p(:,1:Ny-1)); 

end
