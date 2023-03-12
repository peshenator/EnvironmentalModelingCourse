% Add the gradients of P' (Pprime) to the momemeta
function [ustar,vstar] = MomentumPressureGrad(ustar,vstar,Pstar)
global Nx Ny dx dy dt;

dtdx = dt/dx;
dtdy = dt/dy;
% add x-grad of Pstar to the u-velocity
for i=1:Nx+1
    for j=1:Ny
        if ( i==1 )
            ustar(i,j) = ustar(i,j) - dtdx*( Pstar(i,j) - Pstar(Nx,j) );
        elseif( i==Nx+1 )
            ustar(i,j) = ustar(i,j) - dtdx*( Pstar(1,j) - Pstar(i-1,j) );
        else
            ustar(i,j) = ustar(i,j) - dtdx*( Pstar(i,j) - Pstar(i-1,j) );
        end
    end
end

% add y-grad of Pstar to the v-velocity
for i=1:Nx
    for j=1:Ny+1
        if ( j==1 )
            vstar(i,j) = vstar(i,j) - dtdy*( Pstar(i,j) - Pstar(i,Ny) );
        elseif( j==Ny+1 )
            vstar(i,j) = vstar(i,j) - dtdy*( Pstar(i,1) - Pstar(i,j-1) );
        else
            vstar(i,j) = vstar(i,j) - dtdy*( Pstar(i,j) - Pstar(i,j-1) );
        end
    end
end


end
