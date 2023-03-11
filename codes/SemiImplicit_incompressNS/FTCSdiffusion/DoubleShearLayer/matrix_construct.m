imax = 6;
jmax = 6;
Bx = randn(imax+1,jmax);    % some matrix analogoues to H in your projects, (imax+1)*jmax !
By = randn(imax,jmax+1);    % some matrix analogoues to H in your projects, imax*(jmax+1) !
rho= randn(imax,jmax);
vel= randn(imax,jmax);
dtdx2 = 0.123;
dtdy2 = 0.123;
M = zeros(imax*jmax);       % Matrix that we need to contruct 
RHS = zeros(imax*jmax,1);
i=0;    
j=0;
p=1;                        % p=1 for scalar or p=3 for exaple if eta is a vector of length 3
for k = 1:imax*jmax
    if(rem(k-1,imax)/k==0)
        j=j+1;
        i=1;
    else
        i=i+1;
    end
    r=(p*(k-1)+1:p*k);
	M(r,r) = rho(i,j)*eye(p) - dtdx2*(Bx(i+1,j) + Bx(i,j)) - dtdy2*(By(i,j+1) + By(i,j));
    if (j == 1)
        if (i == 1)
            M(r,r+p)      = dtdx2*Bx(i+1,j);
            M(r,p*imax+r) = dtdy2*By(i,j+1);
        elseif(i == imax)
            M(r,r-p)      = dtdx2*Bx(i,j  );
            M(r,p*imax+r) = dtdy2*By(i,j+1);
        else
            M(r,r-p)      = dtdx2*Bx(i,j  );
            M(r,r+p)      = dtdx2*Bx(i+1,j);
            M(r,p*imax+r) = dtdy2*By(i,j+1);
        end
    elseif(j==jmax)
        if (i == 1)
            M(r,r+p)       = dtdx2*Bx(i+1,j);
            M(r,-p*imax+r) = dtdy2*By(i,j);
        elseif(i == imax)
            M(r,r-p)       = dtdx2*Bx(i,j);
            M(r,-p*imax+r) = dtdy2*By(i,j);
        else
            M(r,r-p)       = dtdx2*Bx(i  ,j);
            M(r,r+p)       = dtdx2*Bx(i+1,j);
            M(r,-p*imax+r) = dtdy2*By(i ,j);
        end
    else
        if(i == 1)
            M(r,r+p)       = dtdx2*Bx(i+1,j);
            M(r,p*imax+r)  = dtdy2*By(i,j+1);
            M(r,-p*imax+r) = dtdy2*By(i ,j);
        elseif(i==imax)
            M(r,r-p)       = dtdx2*Bx(i,j  );
            M(r,p*imax+r)  = dtdy2*By(i,j+1);
            M(r,-p*imax+r) = dtdy2*By(i ,j);
        else
            M(r,r-p)       = dtdx2*Bx(i  ,j);
            M(r,r+p)       = dtdx2*Bx(i+1,j);
            M(r,-p*imax+r) = dtdy2*By(i,j  );
            M(r,p*imax+r)  = dtdy2*By(i,j+1);
        end
    end
    RHS(r) = vel(i,j);
end