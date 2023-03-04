% Add the gradients of P' (Pprime) to the momemeta
function [u,v] = DivFreeVelocityCorrection(u,v,ustar,vstar,Pprime)
global imax jmax dx dy dt;

dtdx = dt/dx;
dtdy = dt/dy;
% add x-grad of Pstar to the u-velocity
for i=1:imax+1
    for j=1:jmax
        if ( i==1 )
            u(i,j) = 0;
        elseif( i==imax+1 )
            u(i,j) = 0;
        else
            u(i,j) = ustar(i,j) - dtdx*( Pprime(i,j) - Pprime(i-1,j) );
        end
    end
end

% add y-grad of Pstar to the v-velocity
for i=1:imax
    for j=1:jmax+1
        if ( j==1 )
            v(i,j) = 0;
        elseif( j==jmax+1 )
            v(i,j) = 0;
        else
            v(i,j) = vstar(i,j) - dtdy*( Pprime(i,j) - Pprime(i,j-1) );
        end
    end
end


end

