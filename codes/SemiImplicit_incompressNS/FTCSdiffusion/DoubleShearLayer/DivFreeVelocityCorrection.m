% Add the gradients of P' (Pprime) to the momemeta
function [u,v] = DivFreeVelocityCorrection(u,v,ustar,vstar,Pprime)
global Nx Ny dx dy dt;

dtdx = dt/dx;
dtdy = dt/dy;
% add x-grad of Pstar to the u-velocity
for i=1:Nx+1
    for j=1:Ny
        if ( i==1 )
            u(i,j) = ustar(i,j) - dtdx*( Pprime(i,j) - Pprime(Nx,j) );
        elseif( i==Nx+1 )
            u(i,j) = ustar(i,j) - dtdx*( Pprime(1,j) - Pprime(i-1,j) );
        else
            u(i,j) = ustar(i,j) - dtdx*( Pprime(i,j) - Pprime(i-1,j) );
        end
    end
end

% add y-grad of Pstar to the v-velocity
for i=1:Nx
    for j=1:Ny+1
        if ( j==1 )
            v(i,j) = vstar(i,j) - dtdy*( Pprime(i,j) - Pprime(i,Ny) );
        elseif( j==Ny+1 )
            v(i,j) = vstar(i,j) - dtdy*( Pprime(i,1) - Pprime(i,j-1) );
        else
            v(i,j) = vstar(i,j) - dtdy*( Pprime(i,j) - Pprime(i,j-1) );
        end
    end
end


end

