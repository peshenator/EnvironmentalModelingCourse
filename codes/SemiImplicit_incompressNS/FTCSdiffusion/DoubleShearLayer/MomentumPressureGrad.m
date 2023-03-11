% Add the gradients of P' (Pprime) to the momemeta
function [ustar,vstar] = MomentumPressureGrad(ustar,vstar,Pstar)
global imax jmax dx dy dt;

dtdx = dt/dx;
dtdy = dt/dy;
% add x-grad of Pstar to the u-velocity
for i=1:imax+1
    for j=1:jmax
        if ( i==1 )
            ustar(i,j) = ustar(i,j) - dtdx*( Pstar(i,j) - Pstar(imax,j) );
        elseif( i==imax+1 )
            ustar(i,j) = ustar(i,j) - dtdx*( Pstar(1,j) - Pstar(i-1,j) );
        else
            ustar(i,j) = ustar(i,j) - dtdx*( Pstar(i,j) - Pstar(i-1,j) );
        end
    end
end

% add y-grad of Pstar to the v-velocity
for i=1:imax
    for j=1:jmax+1
        if ( j==1 )
            vstar(i,j) = vstar(i,j) - dtdy*( Pstar(i,j) - Pstar(i,jmax) );
        elseif( j==jmax+1 )
            vstar(i,j) = vstar(i,j) - dtdy*( Pstar(i,1) - Pstar(i,j-1) );
        else
            vstar(i,j) = vstar(i,j) - dtdy*( Pstar(i,j) - Pstar(i,j-1) );
        end
    end
end


end
