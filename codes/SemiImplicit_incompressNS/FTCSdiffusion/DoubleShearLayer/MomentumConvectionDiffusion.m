% update ustar,vstar with the convection and diffusion terms
function [ustar,vstar] = MomentumConvectionDiffusion(u,v,T)
global Nx Ny dt dx dy nu  beta g T0;

ustar = u;
vstar = v;

dtdx = dt/dx;
dtdy = dt/dy;

dtdx2 = dt/dx^2;
dtdy2 = dt/dy^2;

for i=1:Nx+1
    for j=1:Ny
        % dU/dX
        if ( i==1 )
            uac = 0.5*(u(i,j) + u(i+1,j));  % average velocity at the cell-centers
            um  = 0.5*(uac - abs(uac));     % left moving 
            uac = 0.5*(u(i,j) + u(Nx,j));  % average velocity at the cell-centers
            up  = 0.5*(uac + abs(uac));     % right moving wave
            ustar(i,j) = ustar(i,j) - ...
                         dtdx*(um*( u(i+1,j)-u(i,j) ) + up*( u(i,j)-u(Nx,j) ))...
                    + nu*dtdx2*(  ( u(i+1,j)-u(i,j) ) -    ( u(i,j)-u(Nx,j) ) );
        elseif(i==Nx+1)
            uac = 0.5*(u(i,j) + u(2,j));  % average velocity at the cell-centers
            um  = 0.5*(uac - abs(uac));     % left moving wave
            uac = 0.5*(u(i,j) + u(i-1,j));  % average velocity at the cell-centers
            up  = 0.5*(uac + abs(uac));     % right moving wave
            ustar(i,j) = ustar(i,j) - ...
                         dtdx*( um*( u(2,j)-u(i,j) ) + up*( u(i,j) - u(i-1,j) ) )...
                    + nu*dtdx2*(   ( u(2,j)-u(i,j) ) -    ( u(i,j) - u(i-1,j) ) );
        else
            uav = 0.5*( u(i+1,j) + u(i,j) );
            um  = 0.5*( uav - abs(uav) );        % u^- in lecture notes 
            uav = 0.5*( u(i,j) + u(i-1,j) );    % compute average vel. at barycenter
            up  = 0.5*( uav + abs(uav) );       % u^+ in lecture notes 
            ustar(i,j) = ustar(i,j) - ...
                         dtdx*( um*( u(i+1,j)-u(i,j) ) + up*( u(i,j)-u(i-1,j) ) )...
                    + nu*dtdx2*(   ( u(i+1,j)-u(i,j) ) -    ( u(i,j)-u(i-1,j) ) );
        end
        % dU/dY
        if ( j==1 )
            vav = 0.5*( v(bcj_1(i-1,Nx),j+1) + v(bcj(i,Nx+1),j+1) );    % average v at the corners
            vm  = 0.5*(vav - abs(vav));
            vav = 0.5*( v(bcj_1(i-1,Nx),j  ) + v(bcj(i,Nx+1),j  ) );    % average v at the corners
            vp  = 0.5*(vav + abs(vav)); 
            ustar(i,j) = ustar(i,j) - ...
                        dtdy*( vm*( u(i,j+1)-u(i,j) ) + vp*( u(i,j)-u(i,Ny) ) ) ...
                   + nu*dtdy2*(   ( u(i,j+1)-u(i,j) ) -    ( u(i,j)-u(i,Ny) ) );
        elseif( j==Ny )
            vav = 0.5*( v(bcj_1(i-1,Nx),j+1) + v(bcj(i,Nx+1),j+1) );    % average v at the corners
            vm  = 0.5*(vav - abs(vav)); 
            vav = 0.5*( v(bcj_1(i-1,Nx),j  ) + v(bcj(i,Nx+1),j  ) );    % average v at the corners
            vp  = 0.5*(vav + abs(vav)); 
            ustar(i,j) = ustar(i,j) - ...
                         dtdy*( vm*( u(i,1)-u(i,j) ) + vp*( u(i,j)-u(i,j-1) ) ) ...
                    + nu*dtdy2*(   ( u(i,1)-u(i,j) ) -    ( u(i,j)-u(i,j-1) ) );
        else
            vav = 0.5*( v(bcj_1(i-1,Nx),j+1) + v(bcj(i,Nx+1),j+1) );    % average v at the corners
            vm  = 0.5*(vav - abs(vav)); 
            vav = 0.5*( v(bcj_1(i-1,Nx),j  ) + v(bcj(i,Nx+1),j  ) );    % average v at the corners
            vp  = 0.5*(vav + abs(vav));
            ustar(i,j) = ustar(i,j) - ...
                         dtdy*( vm*( u(i,j+1)-u(i,j) ) + vp*( u(i,j)-u(i,j-1) ) ) ...
                    + nu*dtdy2*(   ( u(i,j+1)-u(i,j) ) -    ( u(i,j)-u(i,j-1) ) );
        end
    end
end

%% v-component of the velocity:
for i=1:Nx
    for j=1:Ny+1
        % dV/dX
        if(i==1)
            uav = 0.5*( u(i+1,bcj_1(j-1,Ny)) + u(i+1,bcj(j,Ny+1) ) );
            um  = 0.5*( uav-abs(uav) );
            uav = 0.5*( u(i  ,bcj_1(j-1,Ny)) + u(i  ,bcj(j,Ny+1) ) );
            up  = 0.5*(uav + abs(uav));
            vstar(i,j) = vstar(i,j) - ...
                        dtdx*( um*(v(i+1,j)-v(i,j)) + up*(v(i,j)-v(Nx,j)) ) ...
                   + nu*dtdx2*(   (v(i+1,j)-v(i,j)) -    (v(i,j)-v(Nx,j)) );
        elseif(i==Nx)
            uav = 0.5*( u(i+1,bcj_1(j-1,Ny)) + u(i+1,bcj(j,Ny+1) ) );
            um  = 0.5*(uav - abs(uav));
            uav = 0.5*( u(i  ,bcj_1(j-1,Ny)) + u(i  ,bcj(j,Ny+1) ) );
            up  = 0.5*(uav + abs(uav));           
            vstar(i,j) = vstar(i,j) - ...
                        dtdx*( um*(v(1,j)-v(i,j)) + up*(v(i,j)-v(i-1,j)) ) ...
                   + nu*dtdx2*(   (v(1,j)-v(i,j)) -    (v(i,j)-v(i-1,j)) );
        else
            uav = 0.5*( u(i+1,bcj_1(j-1,Ny) ) + u(i+1,bcj(j,Ny+1) ) ); 
            um  = 0.5*( uav-abs(uav) );
            uav = 0.5*( u(i  ,bcj_1(j-1,Ny) ) + u(i  ,bcj(j,Ny+1) ) ); 
            up  = 0.5*( uav+abs(uav) );
            if (i > 1 && j > 3)
                um = um;
            end
            vstar(i,j) = vstar(i,j) ...
                         - dtdx*( um*(v(i+1,j)-v(i,j)) + up*(v(i,j)-v(i-1,j)) ) ... 
                     + nu*dtdx2*(    (v(i+1,j)-v(i,j)) -    (v(i,j)-v(i-1,j)) ); 
        end
        % dV/dY
        if(j==1)
            vac = 0.5*(v(i,j+1)+v(i,j));
            vm  = 0.5*(vac - abs(vac));
            vac = 0.5*(v(i,j) + v(i,Ny));
            vp  = 0.5*(vac + abs(vac));
            vstar(i,j) = vstar(i,j) - ...
                         dtdy*( vm*(v(i,j+1)-v(i,j)) + vp*(v(i,j)-v(i,Ny)) ) ... 
                    + nu*dtdy2*(   (v(i,j+1)-v(i,j)) -    (v(i,j)-v(i,Ny)) );
        elseif(j==Ny+1)
            vac = 0.5*(v(i,2)+v(i,j));
            vm  = 0.5*(vac - abs(vac));
            vac = 0.5*(v(i,j) + v(i,j-1));
            vp  = 0.5*(vac + abs(vac)); 
            vstar(i,j) = vstar(i,j) - ...
                        dtdy*( vm*(v(i,2)-v(i,j)) + vp*(v(i,j  )-v(i,j-1)) ) ... 
                   + nu*dtdy2*(   (v(i,2)-v(i,j)) -    (v(i,j  )-v(i,j-1)) );
        else
            vac=0.5*(v(i,j+1)+v(i,j)); 
            vm=0.5*(vac-abs(vac)); 
            vac=0.5*(v(i,j)+v(i,j-1)); 
            vp=0.5*(vac+abs(vac)); 
            vstar(i,j) = vstar(i,j) ...
                        - dtdy*( vm*(v(i,j+1)-v(i,j)) + vp*(v(i,j  )-v(i,j-1)) ) ... 
                    + nu*dtdy2*(    (v(i,j+1)-v(i,j)) -    (v(i,j  )-v(i,j-1)) );
            % add the buoyancy forces
%             Tav = 0.5*( T(i,j) + T(i,j-1) ); % average temperature 
%             vstar(i,j) = vstar(i,j) + dt*beta*g*(Tav - T0); 
        end
    end
end




end

