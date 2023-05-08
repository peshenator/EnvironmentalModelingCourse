% update ustar,vstar with the convection and diffusion terms
function [ustar,vstar] = MomConvectionDiffusion2(u,v,T)
global Nx Ny dt dx dy nu uLid beta g T0;

ustar = u;
vstar = v;

dtdx = dt/dx;
dtdy = dt/dy;

dtdx2 = dt/dx^2;
dtdy2 = dt/dy^2;

for i=1:Nx+1
    for j=1:Ny
        % dU/dY
        if ( j==1 )
            vc = 0.5*( v(bcj_1(i-1,Nx),j+1) + v(bcj(i,Nx+1),j+1) );    % average v at the corners
            vm  = 0.5*(vc - abs(vc));
            vp = 0;
            ustar(i,j) = ustar(i,j) - ...
                        dtdy*( vm*( u(i,j+1)-u(i,j) ) + vp*( u(i,j)-0 ) ) ...
                   + nu*dtdy2*(   ( u(i,j+1)-u(i,j) ) -  2*( u(i,j)-0 ) );    % 2 is from  /(dy/2)
        elseif( j==Ny )
            vm  = 0; 
            vc = 0.5*( v(bcj_1(i-1,Nx),j  ) + v(bcj(i,Nx+1),j  ) );    % average v at the corners
            vp  = 0.5*(vc + abs(vc)); 
            ustar(i,j) = ustar(i,j) - ...
                         dtdy*( vm*( uLid-u(i,j) ) + vp*( u(i,j)-u(i,j-1) ) ) ...
                    + nu*dtdy2*( 2*( uLid  -u(i,j) ) -    ( u(i,j)-u(i,j-1) ) );    % 2 is from  /(dy/2)
        else
            vc = 0.5*( v(bcj_1(i-1,Nx),j+1) + v(bcj(i,Nx+1),j+1) );    % average v at the corners
            vm  = 0.5*(vc - abs(vc)); 
            vc = 0.5*( v(bcj_1(i-1,Nx),j  ) + v(bcj(i,Nx+1),j  ) );    % average v at the corners
            vp  = 0.5*(vc + abs(vc));
            ustar(i,j) = ustar(i,j) - ...
                         dtdy*( vm*( u(i,j+1)-u(i,j) ) + vp*( u(i,j)-u(i,j-1) ) ) ...
                    + nu*dtdy2*(   ( u(i,j+1)-u(i,j) ) -    ( u(i,j)-u(i,j-1) ) );
        end
        % dU/dX
        if ( i==1 )
            ustar(i,j) = 0; % wall (no-slip)
        elseif(i==Nx+1)
            ustar(i,j) = 0; % wall (no-slip)
        else
            ub = 0.5*( u(i+1,j) + u(i,j) );     % compute average vel. at barycenter
            um  = 0.5*( ub - abs(ub) );         % u^- in lecture notes 
            ub = 0.5*( u(i,j) + u(i-1,j) );     % compute average vel. at barycenter
            up  = 0.5*( ub + abs(ub) );         % u^+ in lecture notes 
            ustar(i,j) = ustar(i,j) - ...
                         dtdx*( um*( u(i+1,j)-u(i,j) ) + up*( u(i,j)-u(i-1,j) ) )...
                    + nu*dtdx2*(   ( u(i+1,j)-u(i,j) ) -    ( u(i,j)-u(i-1,j) ) );
        end
    end
end

%% v-component of the velocity:
for i=1:Nx
    for j=1:Ny+1
        % dV/dX
        if(i==1)
            uc = 0.5*( u(i+1,bcj_1(j-1,Ny)) + u(i+1,bcj(j,Ny+1) ) );
            um  = 0.5*( uc-abs(uc) );
            up  = 0;
            vstar(i,j) = vstar(i,j) - ...
                        dtdx*( um*(v(i+1,j)-v(i,j)) + up*(v(i,j)-v(Nx,j)) ) ...
                   + nu*dtdx2*(   (v(i+1,j)-v(i,j)) -  2*(v(i,j)-0      ) );
        elseif(i==Nx)
            um  = 0;
            uc = 0.5*( u(i  ,bcj_1(j-1,Ny)) + u(i  ,bcj(j,Ny+1) ) );
            up  = 0.5*(uc + abs(uc));           
            vstar(i,j) = vstar(i,j) - ...
                        dtdx*( um*(v(1,j)-v(i,j)) + up*(v(i,j)-v(i-1,j)) ) ...
                   + nu*dtdx2*( 2*(0     -v(i,j)) -    (v(i,j)-v(i-1,j)) );
        else
            uc = 0.5*( u(i+1,bcj_1(j-1,Ny) ) + u(i+1,bcj(j,Ny+1) ) );   % u at the corners
            um  = 0.5*( uc-abs(uc) );
            uc = 0.5*( u(i  ,bcj_1(j-1,Ny) ) + u(i  ,bcj(j,Ny+1) ) );   % u at the corners
            up  = 0.5*( uc+abs(uc) );
            vstar(i,j) = vstar(i,j) ...
                         - dtdx*( um*(v(i+1,j)-v(i,j)) + up*(v(i,j)-v(i-1,j)) ) ... 
                     + nu*dtdx2*(    (v(i+1,j)-v(i,j)) -    (v(i,j)-v(i-1,j)) ); 
        end
        % dV/dY
        if(j==1)
            vstar(i,j) = 0;
        elseif(j==Ny+1)
            vstar(i,j) = 0;
        else
            vac=0.5*(v(i,j+1)+v(i,j)); 
            vm=0.5*(vac-abs(vac)); 
            vac=0.5*(v(i,j)+v(i,j-1)); 
            vp=0.5*(vac+abs(vac)); 
            vstar(i,j) = vstar(i,j) ...
                        - dtdy*( vm*(v(i,j+1)-v(i,j)) + vp*(v(i,j  )-v(i,j-1)) ) ... 
                    + nu*dtdy2*(    (v(i,j+1)-v(i,j)) -    (v(i,j  )-v(i,j-1)) );
            % add the buoyancy forces
            Tav = 0.5*( T(i,j) + T(i,j-1) ); % average temperature 
            vstar(i,j) = vstar(i,j) + dt*beta*g*(Tav - T0); 
        end
    end
end

%% Now, we add the explciit part of the pressure gradient: 
% [ustar,vstar] = AddGradPressure(ustar,vstar,P,dt);
end

