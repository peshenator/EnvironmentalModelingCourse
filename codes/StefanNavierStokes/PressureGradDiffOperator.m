% Compute the grad(P') and add it to the convection differential operator
function [ustar,vstar] = PressureGradDiffOperator(ustar,vstar,Pstar)
global Nx Ny dx dy dt;

dtdx = dt/dx;
dtdy = dt/dy;
% add x-grad of Pstar to the u-velocity
for i=1:Nx+1
    for j=1:Ny
        if ( i==1 )
            ustar(i,j) = 0;
        elseif( i==Nx+1 )
            ustar(i,j) = 0;
        else
            ustar(i,j) = ustar(i,j) + dtdx*( Pstar(i,j) - Pstar(i-1,j) );
        end
    end
end

% add y-grad of Pstar to the v-velocity
for i=1:Nx
    for j=1:Ny+1
        if ( j==1 )
            vstar(i,j) = 0;
        elseif( j==Ny+1 )
            vstar(i,j) = 0;
        else
            vstar(i,j) = vstar(i,j) + dtdy*( Pstar(i,j) - Pstar(i,j-1) );
        end
    end
end


end
