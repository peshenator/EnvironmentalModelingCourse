% implement the nonlnear advection and diffusion + buoyncy force (Bousinesque approximation)

function [ustar,vstar] = MomentumConvectionDiffusion(u,v,T)
global imax jmax dt dx dy nu uLid beta T0 uWall vWall g;

ustar = u;
vstar = v;

% u-component of the valocity

%SPACE LOOP
for i=1:imax+1
    for j=1:jmax

        % x flux
        if( i == 1)
            ustar(i,j) = 0;         % wall boundary condition
        elseif( i == imax+1)
            ustar(i,j) = 0;         % wall boundary conditon
        else
            uav = 0.5*(u(i,j) + u(i-1,j));  % average velocity at the cell-centers
            up  = 0.5*(uav + abs(uav));     % right moving wave
            uav = 0.5*(u(i,j) + u(i+1,j));  % average velocity at the cell-centers
            um  = 0.5*(uav - abs(uav));     % left moving wave
            ustar(i,j) = ustar(i,j) - dt/dx*(um*(u(i+1,j)-u(i,j)) + up*(u(i,j) - u(i-1,j)));  % convective terms
                                    %+ nu*dt/dx*( (u(i+1,j)-u(i,j))/dx - (u(i,j)-u(i-1,j))/dx ); % diffusive terms
        end

        % y flux
        if(j == 1)
            vav = 0.5*( v(max(1,i-1),j+1) + v(min(imax,i),j+1)); % velocity averages stored at the cell corners
            vm  = 0.5*(vav - abs(vav));     % left going wave
            ustar(i,j) = ustar(i,j) - dt/dy*( vm*(u(i,j+1)-u(i,j)) + vWall )... %; % convective terms
                                    + nu*dt/dy*( 0*(u(i,j+1)-u(i,j))/dy -(0*u(i,j)-uWall)/(dy/2) );            
        elseif(j == jmax)
            vav = 0.5*( v(max(1,i-1),j) + v(min(imax,i),j)); % velocity averages stored at the cell corners
            vp  = 0.5*(vav + abs(vav));     % right going wave
            ustar(i,j) = ustar(i,j) - dt/dy*( vWall + vp*(u(i,j)-u(i,j-1)) )... %; % convective terms
                                    + nu*dt/dy*( (uLid-0*u(i,j))/(dy/2) -0*(u(i,j)-u(i,j-1))/dy );           
        else
            vav = 0.5*( v(max(1,i-1),j+1) + v(min(imax,i),j+1)); % velocity averages stored at the cell corners
            vm  = 0.5*(vav - abs(vav));     % left going wave
            vav = 0.5*( v(max(1,i-1),j) + v(min(imax,i),j)); % velocity averages stored at the cell corners
            vp  = 0.5*(vav + abs(vav));     % right going wave
            ustar(i,j) = ustar(i,j) - dt/dy*( vm*(u(i,j+1)-u(i,j)) + vp*(u(i,j)-u(i,j-1)) ); % convective terms
                                   % + nu*dt/dy*( (u(i,j+1)-u(i,j))/dy -(u(i,j)-u(i,j-1))/dy );
        end
    end
end


% v-component of the valocity
%SPACE LOOP

for i=1:imax
    for j=1:jmax+1
        % x flux
        if(i == 1)
            uav = 0.5*( u(i+1,max(1,j-1)) + u(i+1,min(jmax,j)) );
            um = 0.5*(uav - abs(uav));
            vstar(i,j) = vstar(i,j) - dt/dx*( uWall + um*(v(i+1,j)-v(i,j)) )... %;
                                    + nu*dt/dx*( (v(i+1,j)-v(i,j))/dx - (v(i,j)-vWall)/(dx/2) );
        elseif(i == imax)
            uav = 0.5*( u(i,max(1,j-1)) + u(i,min(jmax,j)) );
            up = 0.5*(uav + abs(uav));
            vstar(i,j) = vstar(i,j) - dt/dx*( up*(v(i,j)-v(i-1,j)) + uWall )... %;
                                    + nu*dt/dx*( (vWall-v(i,j))/(dx/2) - (v(i,j)-v(i-1,j))/dx );
        else
            uav = 0.5*( u(i,max(1,j-1)) + u(i,min(jmax,j)) );
            up = 0.5*(uav + abs(uav));
            uav = 0.5*( u(i+1,max(1,j-1)) + u(i+1,min(jmax,j)) );
            um = 0.5*(uav - abs(uav));
            vstar(i,j) = vstar(i,j) - dt/dx*( up*(v(i,j)-v(i-1,j)) + um*(v(i+1,j)-v(i,j)) );
                                   % + nu*dt/dx*( (v(i+1,j)-v(i,j))/dx - (v(i,j)-v(i-1,j))/dx );
        end
        
        % y flux
        if(j==1)
            vstar(i,j) = 0;     % wall boundary condition
        elseif(j==jmax+1)
            vstar(i,j) = 0;
        else
            vav = 0.5*(v(i,j+1)+v(i,j));
            vm  = 0.5*(vav - abs(vav));
            vav = 0.5*(v(i,j) + v(i,j-1));
            vp  = 0.5*(vav + abs(vav));
            vstar(i,j) = vstar(i,j) - dt/dy*( vp*(v(i,j)-v(i,j-1)) + vm*(v(i,j+1)-v(i,j)) );
                                   % +nu*dt/dy*( (v(i,j+1)-v(i,j))/dy - (v(i,j)-v(i,j-1))/dy );
            
            % add the buyoncy forces to the v velocity
            Tav = 0.5*(T(i,j)+T(i,j-1)); % T cell centre -> average
            vstar(i,j) = vstar(i,j) + dt*beta*g*(Tav - T0);
        end
    end
end

end