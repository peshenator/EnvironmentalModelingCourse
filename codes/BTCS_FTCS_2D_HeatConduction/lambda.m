% define thermal diffusivity lambda = K/(rho*C) at the cell centers
% and the cell faces
function [L,Lxm,Lxp,Lym,Lyp] = lambda(L,Lxm,Lxp,Lym,Lyp,x,y)
global Nx Ny;
    
    % heat diffusifity coefficient in the cell centers
    for i = 1:length(x)
        for j = 1:length(y)
            if (x(i) < (x(1) + x(end))/2)
            % if ((x(i) - 0)^2 +(y(j) - 0)^2 < 0.5)
                L(i,j) = 1;
            else
                L(i,j) = 5;
            end
        end
    end

    Lxm(1     ,:) = L(1,:);
    Lxm(2:Nx  ,:) = 0.5*( L(2:Nx,:) + L(1:Nx-1,:) );
    Lxp(1:Nx-1,:) = Lxm(2:Nx,:);
    Lxp(Nx    ,:) = L(Nx,:);
    
    Lym(:,1     ) = L(:,1);
    Lym(:,2:Ny  ) = 0.5*( L(:,2:Ny) + L(:,1:Ny-1) );
    Lyp(:,1:Ny-1) = Lym(:,2:Ny);
    Lyp(:,Ny    ) = L(:,Ny);
end