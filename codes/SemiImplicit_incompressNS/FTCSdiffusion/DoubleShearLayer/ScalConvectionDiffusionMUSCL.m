% Convection + diffusion of the temperature field

%         i+1/2
%    |     |     |
%    |   Tm|Tp   |
%    |     |     |
%    |--x--|--x--|
%      i     i+1 
function Tnew = ScalConvectionDiffusionMUSCL(u,v,T,dt)

global Nx Ny dx dy lambda;

Tnew = T;

dtdx = dt/dx;
dtdy = dt/dy;


%% Piecewise linear reconstruction in space and time
Txp = zeros(Nx,Ny);
Txm = zeros(Nx,Ny);
Typ = zeros(Nx,Ny);
Tym = zeros(Nx,Ny);
slopeX = zeros(Nx,Ny);
slopeY = zeros(Nx,Ny);
for i=1:Nx
    for j=1:Ny
        % X - direction
        if (i==1)
            slopeX(i,j) = minmod(T(i,j)-T(Nx,j),T(i+1,j)-T(i,j));
        elseif(i==Nx)
            slopeX(i,j) = minmod(T(i,j)-T(i-1,j),T(1,j)-T(i,j));
        else
            slopeX(i,j) = minmod(T(i,j)-T(i-1,j),T(i+1,j)-T(i,j));
        end
        %
        % Y - direction
        if (j==1)
            slopeY(i,j) = minmod(T(i,j)-T(i,Ny),T(i,j+1)-T(i,j));
        elseif(j==Ny)
            slopeY(i,j) = minmod(T(i,j)-T(i,j-1),T(i,1)-T(i,j));
        else
            slopeY(i,j) = minmod(T(i,j)-T(i,j-1),T(i,j+1)-T(i,j));
        end
        % Compute extrapolated X-data at the cell boundaries
        Txp(i,j) = T(i,j) - 0.5*slopeX(i,j);
        Txm(i,j) = T(i,j) + 0.5*slopeX(i,j);
        % Compute extrapolated Y-data at the cell boundaries
        Typ(i,j) = T(i,j) - 0.5*slopeY(i,j);
        Tym(i,j) = T(i,j) + 0.5*slopeY(i,j);
        % Compute time derivative from the PDE:
        T_t =-( u(i+1,j)*Txm(i,j) - u(i,j)*Txp(i,j) )/dx - ( v(i,j+1)*Tym(i,j) - v(i,j)*Typ(i,j) )/dy;
        % Evolve in time boundary extrapolated X-data
        Txp(i,j) = Txp(i,j) + 0.5*dt*T_t;
        Txm(i,j) = Txm(i,j) + 0.5*dt*T_t;
        % Evolve in time boundary extrapolated Y-data
        Typ(i,j) = Typ(i,j) + 0.5*dt*T_t;
        Tym(i,j) = Tym(i,j) + 0.5*dt*T_t;
    end
end


for i=1:Nx
    for j=1:Ny
        % x-flux
        if ( i==1 )
            % advection:
            fm = 0.5*u(i  ,j)*( Txp(i  ,j) + Txm(Nx,j) ) - 0.5*abs(u(i  ,j))*( Txp(i  ,j) - Txm(Nx,j) );
            fp = 0.5*u(i+1,j)*( Txp(i+1,j) + Txm(i ,j) ) - 0.5*abs(u(i+1,j))*( Txp(i+1,j) - Txm(i ,j) );
             % diffusion lambda*T_x (Fourier law):
%             fm = fm - lambda*( T(i  ,j) - T(Nx,j) )/dx;
%             fp = fp - lambda*( T(i+1,j) - T(i ,j) )/dx;
        elseif( i==Nx )
            fm = 0.5*u(i   ,j)*( Txp(i,j) + Txm(i-1,j) ) - 0.5*abs(u(i,j))*( Txp(i,j) - Txm(i-1,j) );
            fp = 0.5*u(Nx+1,j)*( Txp(1,j) + Txm(i  ,j) ) - 0.5*abs(u(1,j))*( Txp(1,j) - Txm(i  ,j) );
%             fm = fm - lambda*( T(i,j) - T(i-1,j) )/dx;
%             fp = fp - lambda*( T(1,j) - T(i  ,j) )/dx;
        else
            fm = 0.5*u(i  ,j)*( Txp(i  ,j) + Txm(i-1,j) ) - 0.5*abs(u(i  ,j))*( Txp(i  ,j) - Txm(i-1,j) );
            fp = 0.5*u(i+1,j)*( Txp(i+1,j) + Txm(i  ,j) ) - 0.5*abs(u(i+1,j))*( Txp(i+1,j) - Txm(i  ,j) );
%             fm = fm - lambda*( T(i  ,j) - T(i-1,j) )/dx;
%             fp = fp - lambda*( T(i+1,j) - T(i  ,j) )/dx;
        end
       % y-fluxes 
        if(j==1)
            gm = 0.5*v(i,j  )*( Typ(i,j  ) + Tym(i,Ny) ) - 0.5*abs(v(i,j  ))*(Typ(i,j  ) - Tym(i,Ny) ); 
            gp = 0.5*v(i,j+1)*( Typ(i,j+1) + Tym(i,j ) ) - 0.5*abs(v(i,j+1))*(Typ(i,j+1) - Tym(i,j ) );
%             gm = gm - lambda*( T(i,j  ) - T(i,Ny) )/dy; 
%             gp = gp - lambda*( T(i,j+1) - T(i,j ) )/dy;
        elseif(j==Ny)
            gm = 0.5*v(i,j)*( Typ(i,j) + Tym(i,j-1) ) - 0.5*abs(v(i,j))*( Typ(i,j) - Tym(i,j-1) ); 
            gp = 0.5*v(i,1)*( Typ(i,1) + Tym(i,j  ) ) - 0.5*abs(v(i,1))*( Typ(i,1) - Tym(i,j  ) );
%             gm = gm - lambda*( T(i,j) - T(i,j-1) )/dy; 
%             gp = gp - lambda*( T(i,1) - T(i,j  ) )/dy;
        else
            gm = 0.5*v(i,j  )*( Typ(i,  j) + Tym(i,j-1) ) - 0.5*abs(v(i,j  ))*( Typ(i,j  ) - Tym(i,j-1) ); 
            gp = 0.5*v(i,j+1)*( Typ(i,j+1) + Tym(i,j  ) ) - 0.5*abs(v(i,j+1))*( Typ(i,j+1) - Tym(i,j  ) );
%             gm = gm - lambda*( T(i,j  ) - T(i,j-1) )/dy; 
%             gp = gp - lambda*( T(i,j+1) - T(i,j  ) )/dy;
        end
        Tnew(i,j) = T(i,j) - dtdx*( fp - fm ) - dtdy*( gp - gm );
    end
end


end