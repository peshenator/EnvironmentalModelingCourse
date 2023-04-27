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
        % y flux, dU/dY
        if(j == 1)
            vav = 0.5*( v(bcj_1(i-1,Nx),j+1) + v(bcj(i,Nx+1),j+1) );    % average v at the corners
            vm  = 0.5*(vav - abs(vav));
            vav = 0.5*( v(bcj_1(i-1,Nx),j  ) + v(bcj(i,Nx+1),j  ) );    % average v at the corners
            vp  = 0.5*(vav + abs(vav)); 
            ustar(i,j) = ustar(i,j) - dtdy*( vm*( u     (i,j+1)-u(i,j) ) + vp*( u(i,j)-u(i,Ny) ) );
        elseif( j==Ny )
            vav = 0.5*( v(bcj_1(i-1,Nx),j+1) + v(bcj(i,Nx+1),j+1) );    % average v at the corners
            vm  = 0.5*(vav - abs(vav)); 
            vav = 0.5*( v(bcj_1(i-1,Nx),j  ) + v(bcj(i,Nx+1),j  ) );    % average v at the corners
            vp  = 0.5*(vav + abs(vav)); 
            ustar(i,j) = ustar(i,j) - dtdy*( vm*( u(i,1)-u(i,j) ) + vp*( u(i,j)-u(i,j-1) ) );
        else
            vav = 0.5*( v(bcj_1(i-1,Nx),j+1) + v(bcj(i,Nx+1),j+1) );    % average v at the corners
            vm  = 0.5*(vav - abs(vav)); 
            vav = 0.5*( v(bcj_1(i-1,Nx),j  ) + v(bcj(i,Nx+1),j  ) );    % average v at the corners
            vp  = 0.5*(vav + abs(vav));
            ustar(i,j) = ustar(i,j) - dtdy*( vm*( u(i,j+1)-u(i,j) ) + vp*( u(i,j)-u(i,j-1) ) );
        end
        % x flux
        if( i == 1)
            uac = 0.5*(u(i,j) + u(i+1,j));  % average velocity at the cell-centers
            um  = 0.5*(uac - abs(uac));     % left moving 
            uac = 0.5*(u(i,j) + u(Nx,j));  % average velocity at the cell-centers
            up  = 0.5*(uac + abs(uac));     % right moving wave
            ustar(i,j) = ustar(i,j) - dtdx*(um*( u(i+1,j)-u(i,j) ) + up*( u(i,j)-u(Nx,j) ));
        elseif(i==Nx+1)
            uac = 0.5*(u(i,j) + u(2,j));  % average velocity at the cell-centers
            um  = 0.5*(uac - abs(uac));     % left moving wave
            uac = 0.5*(u(i,j) + u(i-1,j));  % average velocity at the cell-centers
            up  = 0.5*(uac + abs(uac));     % right moving wave
            ustar(i,j) = ustar(i,j) - dtdx*( um*( u(2,j)-u(i,j) ) + up*( u(i,j) - u(i-1,j) ) );
        else
            uav = 0.5*(u(i,j) + u(i+1,j));  % average velocity at the cell-centers
            um  = 0.5*(uav - abs(uav));     % left moving wave
            uav = 0.5*(u(i,j) + u(i-1,j));  % average velocity at the cell-centers
            up  = 0.5*(uav + abs(uav));     % right moving wave
            ustar(i,j) = ustar(i,j) - dtdx*(um*(u(i+1,j)-u(i,j)) + up*(u(i,j) - u(i-1,j)));  % convective terms
        end
    end
end

%% v-component of the valocity
for i=1:Nx
    for j=1:Ny+1
        % x flux
        if(i == 1)
            uav = 0.5*( u(i+1,bcj_1(j-1,Ny)) + u(i+1,bcj(j,Ny+1) ) );
            um  = 0.5*( uav-abs(uav) );
            uav = 0.5*( u(i  ,bcj_1(j-1,Ny)) + u(i  ,bcj(j,Ny+1) ) );
            up  = 0.5*(uav + abs(uav));
            vstar(i,j) = vstar(i,j) - dtdx*( um*(v(i+1,j)-v(i,j)) + up*(v(i,j)-v(Nx,j)) );
        elseif(i==Nx)
            uav = 0.5*( u(i+1,bcj_1(j-1,Ny)) + u(i+1,bcj(j,Ny+1) ) );
            um  = 0.5*(uav - abs(uav));
            uav = 0.5*( u(i  ,bcj_1(j-1,Ny)) + u(i  ,bcj(j,Ny+1) ) );
            up  = 0.5*(uav + abs(uav));           
            vstar(i,j) = vstar(i,j) - dtdx*( um*(v(1,j)-v(i,j)) + up*(v(i,j)-v(i-1,j)) );
        else
            uav = 0.5*( u(i+1,bcj_1(j-1,Ny) ) + u(i+1,bcj(j,Ny+1) ) ); 
            um  = 0.5*( uav-abs(uav) );
            uav = 0.5*( u(i  ,bcj_1(j-1,Ny) ) + u(i  ,bcj(j,Ny+1) ) ); 
            up  = 0.5*( uav+abs(uav) );
            vstar(i,j) = vstar(i,j) - dtdx*( um*(v(i+1,j)-v(i,j)) + up*(v(i,j)-v(i-1,j)) ); 
        end
        % y flux
        if(j==1)
            vstar(i,j) = 0;     % wall boundary condition
        elseif(j==Ny+1)
            vstar(i,j) = 0;
        else
            vav = 0.5*(v(i,j+1)+v(i,j));
            vm  = 0.5*(vav - abs(vav));
            vav = 0.5*(v(i,j) + v(i,j-1));
            vp  = 0.5*(vav + abs(vav));
            vstar(i,j) = vstar(i,j) - dtdy*( vp*(v(i,j)-v(i,j-1)) + vm*(v(i,j+1)-v(i,j)) );
        end
    end
end

end