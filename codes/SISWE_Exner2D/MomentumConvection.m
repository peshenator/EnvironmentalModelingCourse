% implement the nonlnear advection and diffusion + buoyncy force (Bousinesque approximation)

function [ustar,vstar] = MomentumConvection(u,v)
global Nx Ny dt dx dy uWall;

ustar = u;
vstar = v;

dtdx = dt/dx;
dtdy = dt/dy;

%% u-component of the valocity
for i=1:Nx+1
    for j=1:Ny
        % V*dU/dY
        if(j == 1)
            vc = 0.5*( v(bcj_1(i-1,Nx),j+1) + v(bcj(i,Nx+1),j+1) );    % average v at the corners
            vm  = 0.5*(vc - abs(vc));
            vc = 0.5*( v(bcj_1(i-1,Nx),j  ) + v(bcj(i,Nx+1),j  ) );    % average v at the corners
            vp  = 0.5*(vc + abs(vc)); 
            ustar(i,j) = ustar(i,j) - dtdy*( vm*( u     (i,j+1)-u(i,j) ) + vp*( u(i,j)-u(i,Ny) ) );
        elseif( j==Ny )
            vc = 0.5*( v(bcj_1(i-1,Nx),j+1) + v(bcj(i,Nx+1),j+1) );    % average v at the corners
            vm  = 0.5*(vc - abs(vc)); 
            vc = 0.5*( v(bcj_1(i-1,Nx),j  ) + v(bcj(i,Nx+1),j  ) );    % average v at the corners
            vp  = 0.5*(vc + abs(vc)); 
            ustar(i,j) = ustar(i,j) - dtdy*( vm*( u(i,1)-u(i,j) ) + vp*( u(i,j)-u(i,j-1) ) );
        else
            vc = 0.5*( v(bcj_1(i-1,Nx),j+1) + v(bcj(i,Nx+1),j+1) );    % average v at the corners
            vm  = 0.5*(vc - abs(vc)); 
            vc = 0.5*( v(bcj_1(i-1,Nx),j  ) + v(bcj(i,Nx+1),j  ) );    % average v at the corners
            vp  = 0.5*(vc + abs(vc));
            ustar(i,j) = ustar(i,j) - dtdy*( vm*( u(i,j+1)-u(i,j) ) + vp*( u(i,j)-u(i,j-1) ) );
        end
        % U*dU/dX
        if( i == 1)
            ub = 0.5*(u(i,j) + u(i+1,j));  % average velocity at the cell-centers
            um  = 0.5*(ub - abs(ub));     % left moving 
            ub = 0.5*(u(i,j) + u(Nx,j));  % average velocity at the cell-centers
            up  = 0.5*(ub + abs(ub));     % right moving wave
            ustar(i,j) = ustar(i,j) - dtdx*(um*( u(i+1,j)-u(i,j) ) + up*( u(i,j)-u(Nx,j) ));
        elseif(i==Nx+1)
            ub = 0.5*(u(i,j) + u(2,j));  % average velocity at the cell-centers
            um  = 0.5*(ub - abs(ub));     % left moving wave
            ub = 0.5*(u(i,j) + u(i-1,j));  % average velocity at the cell-centers
            up  = 0.5*(ub + abs(ub));     % right moving wave
            ustar(i,j) = ustar(i,j) - dtdx*( um*( u(2,j)-u(i,j) ) + up*( u(i,j) - u(i-1,j) ) );
        else
            ub = 0.5*(u(i,j) + u(i+1,j));  % average velocity at the cell-centers
            um  = 0.5*(ub - abs(ub));     % left moving wave
            ub = 0.5*(u(i,j) + u(i-1,j));  % average velocity at the cell-centers
            up  = 0.5*(ub + abs(ub));     % right moving wave
            ustar(i,j) = ustar(i,j) - dtdx*(um*(u(i+1,j)-u(i,j)) + up*(u(i,j) - u(i-1,j)));  % convective terms
        end
    end
end

%% v-component of the valocity
for i=1:Nx
    for j=1:Ny+1
        % x flux
        if(i == 1)
            uc = 0.5*( u(i+1,bcj_1(j-1,Ny)) + u(i+1,bcj(j,Ny+1) ) );
            um  = 0.5*( uc-abs(uc) );
            uc = 0.5*( u(i  ,bcj_1(j-1,Ny)) + u(i  ,bcj(j,Ny+1) ) );
            up  = 0.5*(uc + abs(uc));
            vstar(i,j) = vstar(i,j) - dtdx*( um*(v(i+1,j)-v(i,j)) + up*(v(i,j)-v(Nx,j)) );
        elseif(i==Nx)
            uc = 0.5*( u(i+1,bcj_1(j-1,Ny)) + u(i+1,bcj(j,Ny+1) ) );
            um  = 0.5*(uc - abs(uc));
            uc = 0.5*( u(i  ,bcj_1(j-1,Ny)) + u(i  ,bcj(j,Ny+1) ) );
            up  = 0.5*(uc + abs(uc));           
            vstar(i,j) = vstar(i,j) - dtdx*( um*(v(1,j)-v(i,j)) + up*(v(i,j)-v(i-1,j)) );
        else
            uc = 0.5*( u(i+1,bcj_1(j-1,Ny) ) + u(i+1,bcj(j,Ny+1) ) ); 
            um  = 0.5*( uc-abs(uc) );
            uc = 0.5*( u(i  ,bcj_1(j-1,Ny) ) + u(i  ,bcj(j,Ny+1) ) ); 
            up  = 0.5*( uc+abs(uc) );
            vstar(i,j) = vstar(i,j) - dtdx*( um*(v(i+1,j)-v(i,j)) + up*(v(i,j)-v(i-1,j)) ); 
        end
        % y flux
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