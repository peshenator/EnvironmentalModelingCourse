function eqP = PoissonEqPressure(P,Nx,Ny,dx2,dy2)
    Nxy = Nx*Ny;
    eqP = zeros(Nxy,1);

    i = 0;
    j = 0;
    for k = 1:Nxy
        if (rem(k-1,Nx)/k == 0)
            j = j+1;
            i = 1;
        else
            i = i+1;
        end
        
        if     (i-1 == 0)
            im = Nx;
            ip = i+1;
        elseif (i+1 == Nx+1)
            im = i-1;
            ip = 1;
        else
            im = i-1;
            ip = i+1;
        end

        if     (j-1 == 0)
            jm = Ny;
            jp = j+1;
        elseif (j+1 == Ny+1)
            jm = j-1;
            jp = 1;
        else
            jm = j-1;
            jp = j+1;
        end

        eqP(k) = (P(im,j ) - 2*P(i,j) + P(ip,j ))/dx2 + ...
                 (P(i ,jm) - 2*P(i,j) + P(i ,jp))/dy2;
    end

end