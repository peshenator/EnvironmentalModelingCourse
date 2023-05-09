% implement the nonlnear advection and diffusion + buoyncy force (Bousinesque approximation)
function [ustar,vstar] = MomConvectionDiffusion(u,v,T)
global Nx Ny dt dx dy nu uLid beta T0 uWall vWall g;

ustar = u;
vstar = v;

dtdx = dt/dx;
dtdy = dt/dy;
dtdx2 = dt/dx^2;
dtdy2 = dt/dy^2;
%% u-component of the valocity

%SPACE LOOP
for i=1:Nx+1
    for j=1:Ny

        % U*dU/dX
        if( i == 1)
            ustar(i,j) = 0;         % wall boundary condition
        elseif( i == Nx+1)
            ustar(i,j) = 0;         % wall boundary conditon
        else
            ub = 0.5*(u(i,j) + u(i-1,j));  % average velocity at the barycenters
            up = 0.5*(ub + abs(ub));     % right moving wave
            ub = 0.5*(u(i,j) + u(i+1,j));  % average velocity at the barycenters
            um = 0.5*(ub - abs(ub));     % left moving wave
            ustar(i,j) = ustar(i,j) - dtdx*( um*( u(i+1,j) - u(i,j) ) + up*( u(i,j) - u(i-1,j) ) );  % convective terms
                                    %+ nu*dt/dx*( (u(i+1,j)-u(i,j))/dx - (u(i,j)-u(i-1,j))/dx ); % diffusive terms
        end

        % V*dU/dY
        if(j == 1)
            vc = 0.5*( v(max(1,i-1),j+1) + v(min(Nx,i),j+1) ); % velocity averages stored at the cell corners
            vm = 0.5*(vc - abs(vc));     % left going wave
            ustar(i,j) = ustar(i,j) - dtdy*( vm*( u(i,j+1) - u(i,j) ) +                vWall   )... %; % convective terms
                                + nu*dtdy2*(  0*( u(i,j+1) - u(i,j) ) - 2*( 0*u(i,j) - uWall ) );            
        elseif(j == Ny)
            vc = 0.5*( v(max(1,i-1),j) + v(min(Nx,i),j) ); % velocity averages stored at the cell corners
            vp = 0.5*(vc + abs(vc));     % right going wave
            ustar(i,j) = ustar(i,j) - dtdy*(     vWall             + vp*( u(i,j) - u(i,j-1) ) )... % ; % convective terms
                                + nu*dtdy2*( 2*( uLid - 0*u(i,j) ) -  0*( u(i,j) - u(i,j-1) ) );           
        else
            vc = 0.5*( v(max(1,i-1),j+1) + v(min(Nx,i),j+1) ); % velocity averages stored at the cell corners
            vm = 0.5*(vc - abs(vc));     % left going wave
            vc = 0.5*( v(max(1,i-1),j  ) + v(min(Nx,i),j  ) ); % velocity averages stored at the cell corners
            vp = 0.5*(vc + abs(vc));     % right going wave
            ustar(i,j) = ustar(i,j) - dtdy*( vm*( u(i,j+1) - u(i,j) ) + vp*( u(i,j) - u(i,j-1) ) ); % convective terms
                              % + nu*dtdy2*(    ( u(i,j+1) - u(i,j) ) -    ( u(i,j) - u(i,j-1) ) );
        end
    end
end


%% v-component of the valocity
%SPACE LOOP

for i=1:Nx
    for j=1:Ny+1
        % U*dV/dX
        if(i == 1)
            uc = 0.5*( u(i+1,max(1,j-1)) + u(i+1,min(Ny,j)) );
            um = 0.5*(uc - abs(uc));
            vstar(i,j) = vstar(i,j) - dtdx*( uWall                 + um*( v(i+1,j) - v(i,j) ) )... %;
                                + nu*dtdx2*( ( v(i+1,j) - v(i,j) ) -  2*( v(i  ,j) - vWall  ) );
        elseif(i == Nx)
            uc = 0.5*( u(i,max(1,j-1)) + u(i,min(Ny,j)) );
            up = 0.5*(uc + abs(uc));
            vstar(i,j) = vstar(i,j) - dtdx*( up*( v(i,j) - v(i-1,j) ) + uWall                 )... %;
                                + nu*dtdx2*(  2*( vWall  - v(i  ,j) ) - ( v(i,j) - v(i-1,j) ) );
        else
            uc = 0.5*( u(i  ,max(1,j-1)) + u(i,min(Ny,j)) );
            up = 0.5*(uc + abs(uc));
            uc = 0.5*( u(i+1,max(1,j-1)) + u(i+1,min(Ny,j)) );
            um = 0.5*(uc - abs(uc));
            vstar(i,j) = vstar(i,j) - dtdx*( up*( v(i  ,j) - v(i-1,j) ) + um*( v(i+1,j) - v(i  ,j) ) );
                              % + nu*dtdx2*(    ( v(i+1,j) - v(i  ,j) ) -    ( v(i  ,j) - v(i-1,j) ) );
        end
        
        % U*dV/dX
        if(j==1)
            vstar(i,j) = 0;     % wall boundary condition
        elseif(j==Ny+1)
            vstar(i,j) = 0;
        else
            vb = 0.5*( v(i,j+1) + v(i,j  ) );
            vm = 0.5*(vb - abs(vb));
            vb = 0.5*( v(i,j  ) + v(i,j-1) );
            vp = 0.5*(vb + abs(vb));
            vstar(i,j) = vstar(i,j) - dtdy*( vp*( v(i,j  ) - v(i,j-1) ) + vm*( v(i,j+1) - v(i,j  ) ) );
                              % + nu*dtdy2*(    ( v(i,j+1) - v(i,j  ) ) -    ( v(i,j  ) - v(i,j-1) ) );
            
            % add the buyoncy forces to the v velocity
            Tav = 0.5*( T(i,j) + T(i,j-1) ); % T cell centre -> average
            vstar(i,j) = vstar(i,j) + dt*beta*g*(Tav - T0);
        end
    end
end

end