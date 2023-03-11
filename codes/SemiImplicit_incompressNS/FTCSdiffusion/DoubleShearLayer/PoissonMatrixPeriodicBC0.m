function [M,RHS] = PoissonMatrixPeriodicBC(imax,jmax,dx2,dy2,rhs)

Bx = 1;
By = 1;
M = zeros(imax*jmax);       % Matrix that we need to contruct 
dx2 = 1/dx2;
dy2 = 1/dy2;
IJ = imax*jmax;
i=0;    
j=0;
p=1;                        % p=1 for scalar or p=3 for exaple if eta is a vector of length 3
for k = 1:IJ
    if(rem(k-1,imax)/k==0)
        j=j+1;
        i=1;
    else
        i=i+1;
    end
    r=(p*(k-1)+1:p*k);
	M(r,r) = - dx2*(Bx + Bx) - dy2*(By + By);
    if (j == 1)
        if (i == 1)
            M(r,r+p)      = dx2*Bx;
            M(r,p*imax)   = dx2*Bx;
            M(r,p*imax+r) = dy2*By;
            M(r,(jmax-1)*imax+r) = dy2*By;  % periodic 
        elseif(i == imax)
            M(r,r-p)      = dx2*Bx;
            M(r,r+p)      = dx2*Bx;         % periodic
            M(r,p*imax+r) = dy2*By;
            M(r,jmax*imax) = dy2*By;  % periodic
        else
            M(r,r-p)      = dx2*Bx;
            M(r,r+p)      = dx2*Bx;
            M(r,p*imax+r) = dy2*By;
            M(r,(jmax-1)*imax+r) = dy2*By;  % periodic
        end
    elseif(j==jmax)
        if (i == 1)
            M(r,r+p)       = dx2*Bx;
            M(r,jmax*imax) = dx2*Bx;  % periodic 
            M(r,-p*imax+r) = dy2*By;
            M(r,p        ) = dy2*By;  % periodic
        elseif(i == imax)
            M(r,r-p)       = dx2*Bx;
            M(r,r-imax+p) = dx2*Bx;  % periodic 
            M(r,-p*imax+r) = dy2*By;
            M(r,jmax        ) = dy2*By;  % periodic
        else
            M(r,r-p)       = dx2*Bx;
            M(r,r+p)       = dx2*Bx;
            M(r,-p*imax+r) = dy2*By;
%             M(r,jmax-        ) = dy2*By;  % periodic
        end
    else
        if(i == 1)
            M(r,r+p)       = dx2*Bx;
            M(r,p*imax+r)  = dy2*By;
            M(r,-p*imax+r) = dy2*By;
        elseif(i==imax)
            M(r,r-p)       = dx2*Bx;
            M(r,p*imax+r)  = dy2*By;
            M(r,-p*imax+r) = dy2*By;
        else
            M(r,r-p)       = dx2*Bx;
            M(r,r+p)       = dx2*Bx;
            M(r,-p*imax+r) = dy2*By;
            M(r,p*imax+r)  = dy2*By;
        end
    end
    RHS(r) = rhs(i,j);
end

end