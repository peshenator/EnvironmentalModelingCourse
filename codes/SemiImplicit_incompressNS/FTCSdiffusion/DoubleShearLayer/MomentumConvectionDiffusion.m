% update ustar,vstar with the convection and diffusion terms
function [ustar,vstar] = MomentumConvectionDiffusion(u,v,T)
global Nx Ny dt dx dy nu  beta g T0;

ustar = u;
vstar = v;

dtdx = dt/dx;
dtdy = dt/dy;

for i=1:Nx+1
    for j=1:Ny
        % x-flux 
        if ( i==1 )
            uac = 0.5*(u(i,j) + u(i+1,j));  % average velocity at the cell-centers
            um  = 0.5*(uac - abs(uac));     % left moving 
            uac = 0.5*(u(i,j) + u(Nx,j));  % average velocity at the cell-centers
            up  = 0.5*(uac + abs(uac));     % right moving wave
            ustar(i,j) = ustar(i,j) - ...
                         dtdx*(um*( u(i+1,j)-u(i,j) ) + up*( u(i,j)-u(Nx,j) ))...
                    + nu*dtdx*(   ( u(i+1,j)-u(i,j) ) -    ( u(i,j)-u(Nx,j) ) )/dx;
        elseif(i==Nx+1)
            uac = 0.5*(u(i,j) + u(2,j));  % average velocity at the cell-centers
            um  = 0.5*(uac - abs(uac));     % left moving wave
            uac = 0.5*(u(i,j) + u(i-1,j));  % average velocity at the cell-centers
            up  = 0.5*(uac + abs(uac));     % right moving wave
            ustar(i,j) = ustar(i,j) - ...
                         dtdx*( um*( u(2,j)-u(i,j) ) + up*( u(i,j) - u(i-1,j) ) )...
                    + nu*dtdx*(    ( u(2,j)-u(i,j) ) -    ( u(i,j) - u(i-1,j) ) )/dx;
        else
            uav = 0.5*( u(i+1,j) + u(i,j) );
            um  = 0.5*( uav - abs(uav) );        % u^- in lecture notes 
            uav = 0.5*( u(i,j) + u(i-1,j) );    % compute average vel. at barycenter
            up  = 0.5*( uav + abs(uav) );       % u^+ in lecture notes 
            ustar(i,j) = ustar(i,j) - ...
                         dtdx*( um*( u(i+1,j)-u(i,j) ) + up*( u(i,j)-u(i-1,j) ) )...
                    + nu*dtdx*(    ( u(i+1,j)-u(i,j) ) -    ( u(i,j)-u(i-1,j) ) )/dx;
        end
        % y-flux
        if ( j==1 )
            vav = 0.5*( v(bcind(i-1,Nx),j+1) +v(min(Nx,i),j+1) );    % average v at the corners
            vm  = 0.5*(vav - abs(vav));
            vav = 0.5*( v(bcind(i-1,Nx),j  ) +v(min(Nx,i),j  ) );    % average v at the corners
            vp  = 0.5*(vav + abs(vav)); 
            ustar(i,j) = ustar(i,j) - ...
                        dtdy*( vm*( u(i,j+1)-u(i,j) ) + vp*( u(i,j)-u(i,Ny) ) ) ...
                   + nu*dtdy*(    ( u(i,j+1)-u(i,j) ) -    ( u(i,j)-u(i,Ny) ) )/dy;
        elseif( j==Ny )
            vav = 0.5*( v(bcind(i-1,Nx),j+1) +v(min(Nx,i),j+1) );    % average v at the corners
            vm  = 0.5*(vav - abs(vav)); 
            vav = 0.5*( v(bcind(i-1,Nx),j) +v(min(Nx,i),j) );    % average v at the corners
            vp  = 0.5*(vav + abs(vav)); 
            ustar(i,j) = ustar(i,j) - ...
                         dtdy*( vm*( u(i,1)-u(i,j) ) + vp*( u(i,j)-u(i,j-1) ) ) ...
                    + nu*dtdy*(    ( u(i,1)-u(i,j) ) -    ( u(i,j)-u(i,j-1) ) )/dy;            
        else
            vav = 0.5*( v(bcind(i-1,Nx),j+1) +v(min(Nx,i),j+1) );    % average v at the corners
            vm  = 0.5*(vav - abs(vav)); 
            vav = 0.5*( v(bcind(i-1,Nx),j) +v(min(Nx,i),j) );    % average v at the corners
            vp  = 0.5*(vav + abs(vav)); 
            ustar(i,j) = ustar(i,j) - ...
                         dtdy*( vm*( u(i,j+1)-u(i,j) ) + vp*( u(i,j)-u(i,j-1) ) ) ...
                    + nu*dtdy*(    ( u(i,j+1)-u(i,j) ) -    ( u(i,j)-u(i,j-1) ) )/dy;
        end
    end
end

%% v-component of the velocity:
for i=1:Nx
    for j=1:Ny+1
        % x-fluxes 
        if(i==1)
            uav = 0.5*( u(i+1,bcind(j-1,Ny)) + u(i+1,min(Ny,j)) );
            um  = 0.5*( uav-abs(uav) );
            uav = 0.5*( u(i  ,bcind(j-1,Ny)) + u(i  ,min(Ny,j)) );
            up  = 0.5*(uav + abs(uav));
            vstar(i,j) = vstar(i,j) - ...
                        dtdx*( um*(v(i+1,j)-v(i,j)) + up*(v(i,j)-v(Nx,j)) ) ...
                   + nu*dtdx*(    (v(i+1,j)-v(i,j)) -    (v(i,j)-v(Nx,j)) )/dx;
        elseif(i==Nx)
            uav = 0.5*( u(i+1,bcind(j-1,Ny)) + u(i+1,min(Ny,j)) );
            um  = 0.5*(uav - abs(uav));
            uav = 0.5*( u(i  ,bcind(j-1,Ny)) + u(i,min(Ny,j)) );
            up  = 0.5*(uav + abs(uav));           
            vstar(i,j) = vstar(i,j) - ...
                        dtdx*( um*(v(1,j)-v(i,j)) + up*(v(i,j)-v(i-1,j)) ) ...
                   + nu*dtdx*(    (v(1,j)-v(i,j)) -    (v(i,j)-v(i-1,j)) )/dx;
        else
            uav = 0.5*( u(i,bcind(1,j-1))+u(i,min(Ny,j)) ); 
            up  = 0.5*( uav+abs(uav) ); 
            uav = 0.5*( u(i+1,bcind(1,j-1))+u(i+1,min(Ny,j)) ); 
            um  = 0.5*( uav-abs(uav) ); 
            vstar(i,j) = vstar(i,j)- ...
                         dtdx*( um*(v(i+1,j)-v(i,j)) + up*(v(i,j)-v(i-1,j)) ) ... 
                    + nu*dtdx*(    (v(i+1,j)-v(i,j)) -    (v(i,j)-v(i-1,j)) )/dx; 
        end
        % y-fluxes 
        if(j==1)
            vac = 0.5*(v(i,j+1)+v(i,j));
            vm  = 0.5*(vac - abs(vac));
            vac = 0.5*(v(i,j) + v(i,Ny));
            vp  = 0.5*(vac + abs(vac));
            vstar(i,j) = vstar(i,j) - ...
                         dtdy*( vm*(v(i,j+1)-v(i,j)) + vp*(v(i,j)-v(i,Ny)) ) ... 
                    + nu*dtdy*(    (v(i,j+1)-v(i,j)) -    (v(i,j)-v(i,Ny)) )/dy;
        elseif(j==Ny+1)
            vac = 0.5*(v(i,2)+v(i,j));
            vm  = 0.5*(vac - abs(vac));
            vac = 0.5*(v(i,j) + v(i,j-1));
            vp  = 0.5*(vac + abs(vac)); 
            vstar(i,j) = vstar(i,j) - ...
                        dtdy*( vm*(v(i,2)-v(i,j)) + vp*(v(i,j  )-v(i,j-1)) ) ... 
                   + nu*dtdy*(    (v(i,2)-v(i,j)) -    (v(i,j  )-v(i,j-1)) )/dy;
        else
            vac=0.5*(v(i,j+1)+v(i,j)); 
            vm=0.5*(vac-abs(vac)); 
            vac=0.5*(v(i,j)+v(i,j-1)); 
            vp=0.5*(vac+abs(vac)); 
            vstar(i,j) = vstar(i,j) - ...
                        dtdy*( vm*(v(i,j+1)-v(i,j)) + vp*(v(i,j  )-v(i,j-1)) ) ... 
                   + nu*dtdy*(    (v(i,j+1)-v(i,j)) -    (v(i,j  )-v(i,j-1)) )/dy;
            % add the buoyancy forces
            Tav = 0.5*( T(i,j) + T(i,j-1) ); % average temperature 
            vstar(i,j) = vstar(i,j) + dt*beta*g*(Tav - T0); 
        end
    end
end




end

