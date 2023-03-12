% Convection + diffusion of the temperature field
function Tnew = TemperatureConvectionDiffusion(u,v,T)

global Nx Ny dt dx dy lambda;

Tnew = T;

dtdx = dt/dx;
dtdy = dt/dy;

for i=1:Nx
    for j=1:Ny
        % x-flux
        if ( i==1 )
            % advection:
            fm = 0.5*u(i  ,j)*( T(i  ,j) + T(Nx,j) ) - 0.5*abs(u(i  ,j))*( T(i  ,j) - T(Nx,j) );
            fp = 0.5*u(i+1,j)*( T(i+1,j) + T(i ,j) ) - 0.5*abs(u(i+1,j))*( T(i+1,j) - T(i ,j) );
             % diffusion lambda*T_x (Fourier law):
            fm = fm - lambda*( T(i  ,j) - T(Nx,j) )/dx;
            fp = fp - lambda*( T(i+1,j) - T(i ,j) )/dx;
        elseif( i==Nx )
            fm = 0.5*u(i,j)*( T(i,j) + T(i-1,j) ) - 0.5*abs(u(i,j))*( T(i,j) - T(i-1,j) );
            fp = 0.5*u(Nx+1,j)*( T(1,j) + T(i  ,j) ) - 0.5*abs(u(1,j))*( T(1,j) - T(i ,j) );
            fm = fm - lambda*( T(i,j) - T(i-1,j) )/dx;
            fp = fp - lambda*( T(1,j) - T(i  ,j) )/dx;
        else
            fm = 0.5*u(i  ,j)*( T(i  ,j) + T(i-1,j) ) - 0.5*abs(u(i  ,j))*( T(i  ,j) - T(i-1,j) );
            fp = 0.5*u(i+1,j)*( T(i+1,j) + T(i  ,j) ) - 0.5*abs(u(i+1,j))*( T(i+1,j) - T(i  ,j) );
            fm = fm - lambda*( T(i  ,j) - T(i-1,j) )/dx;
            fp = fp - lambda*( T(i+1,j) - T(i  ,j) )/dx;
        end
       % y-fluxes 
        if(j==1)
            gm = 0.5*v(i,j)*( T(i,j)+T(i,Ny)) - 0.5*abs(v(i,j))*(T(i,j)-T(i,Ny)); 
            gp = 0.5*v(i,j+1)*( T(i,j+1) + T(i,j) ) - 0.5*abs(v(i,j+1))*(T(i,j+1)-T(i,j));
            gm = gm-lambda*(T(i,j)-T(i,Ny ))/dy; 
            gp = gp-lambda*(T(i,j+1)-T(i,j))/dy;
        elseif(j==Ny)
            gm = 0.5*v(i,j)*(T(i,j)+T(i,j-1)) - 0.5*abs(v(i,j))*(T(i,j)-T(i,j-1)); 
            gp = 0.5*v(i,1)*( T(i,1) + T(i,j) ) - 0.5*abs(v(i,1))*( T(i,1) - T(i,j) );
            gm = gm-lambda*(T(i,j)-T(i,j-1))/dy; 
            gp = gp-lambda*(T(i,1)-T(i,j  ))/dy;
        else
            gm = 0.5*v(i,j  )*( T(i,  j)+T(i,j-1)) - 0.5*abs(v(i,j))*(T(i,j)-T(i,j-1)); 
            gp = 0.5*v(i,j+1)*( T(i,j+1) + T(i,j) ) - 0.5*abs(v(i,j+1))*( T(i,j+1) - T(i,j) );
            gm = gm-lambda*(T(i,j)-T(i,j-1))/dy; 
            gp = gp-lambda*(T(i,j+1)-T(i,j))/dy;
        end
        Tnew(i,j) = T(i,j) - dtdx*( fp - fm ) - dtdy*( gp - gm );
    end
end


end