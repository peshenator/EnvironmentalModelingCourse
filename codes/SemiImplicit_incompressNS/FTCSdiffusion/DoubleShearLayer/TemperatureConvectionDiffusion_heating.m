% Convection + diffusion of the temperature field
function Tnew = TemperatureConvectionDiffusion_heating(u,v,T)

global imax jmax dt dx dy lambda TB xb;

Tnew = T;

dtdx = dt/dx;
dtdy = dt/dy;

for i=1:imax
    for j=1:jmax
        % x-flux
        if ( i==1 )
            fp = 0.5*u(i+1,j)*( T(i+1,j) + T(i  ,j) ) - 0.5*abs(u(i+1,j))*( T(i+1,j) - T(i  ,j) );
            fm = 0;
            fp = fp - lambda*( T(i+1,j) - T(i  ,j) )/dx;    % diffusion lambda*T_x (Fourier law)
        elseif( i==imax )
            fp = 0;
            fm = 0.5*u(i  ,j)*( T(i  ,j) + T(i-1,j) ) - 0.5*abs(u(i  ,j))*( T(i  ,j) - T(i-1,j) );
            fm = fm - lambda*( T(i  ,j) - T(i-1,j) )/dx;    % diffusion lambda*T_x (Fourier law)
        else
            fp = 0.5*u(i+1,j)*( T(i+1,j) + T(i  ,j) ) - 0.5*abs(u(i+1,j))*( T(i+1,j) - T(i  ,j) );
            fm = 0.5*u(i  ,j)*( T(i  ,j) + T(i-1,j) ) - 0.5*abs(u(i  ,j))*( T(i  ,j) - T(i-1,j) );
            fp = fp - lambda*( T(i+1,j) - T(i  ,j) )/dx;    % diffusion lambda*T_x (Fourier law)
            fm = fm - lambda*( T(i  ,j) - T(i-1,j) )/dx;
        end
       % y-fluxes 
        if(j==1)
            if ( TB > 0 )
                temp_bc = TB*exp(-0.5*( xb(i) - 0.5 )^2/0.1^2 );
                gm = -lambda*( T(i,j) - temp_bc)/(dy/2);     % heating at the bottom
            else
                gm = 0;
            end
            gp = 0.5*v(i,j+1)*( T(i,j+1) + T(i,j) ) - 0.5*abs(v(i,j+1))*(T(i,j+1)-T(i,j));
            gp = gp-lambda*(T(i,j+1)-T(i,j))/dy; % diffusion term T_yy            
        elseif(j==jmax)
            gp = 0;
            gm = 0.5*v(i,j)*(T(i,j)+T(i,j-1)) - 0.5*abs(v(i,j))*(T(i,j)-T(i,j-1)); 
            gm = gm-lambda*(T(i,j)-T(i,j-1))/dy; 
        else
            gp = 0.5*v(i,j+1)*( T(i,j+1) + T(i,j) ) - 0.5*abs(v(i,j+1))*( T(i,j+1) - T(i,j) );
            gm = 0.5*v(i,j  )*( T(i,  j)+T(i,j-1)) - 0.5*abs(v(i,j))*(T(i,j)-T(i,j-1)); 
            gp = gp-lambda*(T(i,j+1)-T(i,j))/dy; % diffusion term T_yy 
            gm = gm-lambda*(T(i,j)-T(i,j-1))/dy; 
        end
        Tnew(i,j) = T(i,j) - dtdx*( fp - fm ) - dtdy*( gp - gm );
    end
end


end