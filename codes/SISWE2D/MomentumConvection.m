% implement the nonlnear advection and diffusion + buoyncy force (Bousinesque approximation)

function [ustar,vstar] = MomentumConvection(u,v)
global Nx Ny dt dx dy;

ustar = u;
vstar = v;

dtdx = dt/dx;
dtdy = dt/dy;

vWall = 0;
uWall = 0;
%% u-component of the valocity
for i=1:Nx+1
    for j=1:Ny
        % V*dU/dY
        if(j == 1)
            vc = 0.5*( v(max(1,i-1),j+1) + v(min(Nx,i),j+1)); % velocity averages stored at the cell corners
            vm  = 0.5*(vc - abs(vc));     % left going wave
            ustar(i,j) = ustar(i,j) - dtdy*( vm*(u(i,j+1)-u(i,j)) + vWall );
        elseif(j == Ny)
            vc = 0.5*( v(max(1,i-1),j) + v(min(Nx,i),j)); % velocity averages stored at the cell corners
            vp  = 0.5*(vc + abs(vc));     % right going wave
            ustar(i,j) = ustar(i,j) - dtdy*( vWall + vp*(u(i,j)-u(i,j-1)) );
        else
            vc = 0.5*( v(max(1,i-1),j+1) + v(min(Nx,i),j+1)); % velocity averages stored at the cell corners
            vm  = 0.5*(vc - abs(vc));     % left going wave
            vc = 0.5*( v(max(1,i-1),j) + v(min(Nx,i),j)); % velocity averages stored at the cell corners
            vp  = 0.5*(vc + abs(vc));     % right going wave
            ustar(i,j) = ustar(i,j) - dtdy*( vm*(u(i,j+1)-u(i,j)) + vp*(u(i,j)-u(i,j-1)) ); % convective terms
            % ustar(i,j) = ustar(i,j) - dtdy*( vp*u(i,j) + vm*u(i,j+1) - vp*u(i,j-1) - vm*u(i,j) ); % convective terms
        end
        % U*dU/dX
        if( i == 1)
            ustar(i,j) = 0;         % wall boundary condition
        elseif( i == Nx+1)
            ustar(i,j) = 0;         % wall boundary conditon
        else
            ub = 0.5*(u(i,j) + u(i-1,j));  % average velocity at the cell-centers
            up  = 0.5*(ub + abs(ub));     % right moving wave
            ub = 0.5*(u(i,j) + u(i+1,j));  % average velocity at the cell-centers
            um  = 0.5*(ub - abs(ub));     % left moving wave
            ustar(i,j) = ustar(i,j) - dtdx*(um*(u(i+1,j)-u(i,j)) + up*(u(i,j) - u(i-1,j)));  % convective terms
            % ustar(i,j) = ustar(i,j) - dtdx*( up*u(i,j) + um*u(i+1,j) - (up*u(i-1,j) + um*u(i,j)) );  % convective 
        end
    end
end

%% v-component of the valocity
for i=1:Nx
    for j=1:Ny+1
        % U*dV/dX
        if(i == 1)
            uc = 0.5*( u(i+1,max(1,j-1)) + u(i+1,min(Ny,j)) );
            um = 0.5*(uc - abs(uc));
            vstar(i,j) = vstar(i,j) - dtdx*( uWall + um*(v(i+1,j)-v(i,j)) );
        elseif(i == Nx)
            uc = 0.5*( u(i,max(1,j-1)) + u(i,min(Ny,j)) );
            up = 0.5*(uc + abs(uc));
            vstar(i,j) = vstar(i,j) - dtdx*( up*(v(i,j)-v(i-1,j)) + uWall );
        else
            uc = 0.5*( u(i,max(1,j-1)) + u(i,min(Ny,j)) );
            up = 0.5*(uc + abs(uc));
            uc = 0.5*( u(i+1,max(1,j-1)) + u(i+1,min(Ny,j)) );
            um = 0.5*(uc - abs(uc));
            vstar(i,j) = vstar(i,j) - dtdx*( up*(v(i,j)-v(i-1,j)) + um*(v(i+1,j)-v(i,j)) );
        end
        % V*dV/dY
        if(j==1)
            vstar(i,j) = 0;     % wall boundary condition
        elseif(j==Ny+1)
            vstar(i,j) = 0;
        else
            vb = 0.5*(v(i,j+1)+v(i,j));
            vm  = 0.5*(vb - abs(vb));
            vb = 0.5*(v(i,j) + v(i,j-1));
            vp  = 0.5*(vb + abs(vb));
            vstar(i,j) = vstar(i,j) - dtdy*( vp*(v(i,j)-v(i,j-1)) + vm*(v(i,j+1)-v(i,j)) );
        end
    end
end

end