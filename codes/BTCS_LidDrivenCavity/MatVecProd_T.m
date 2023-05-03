function AT = MatVecProd_T(T)

global dt dx dy lambda Nx Ny;

AT = T;

% SPACE LOOP
for i=1:Nx
    for j=1:Ny

        % the swipe in the x direction
        if (i == 1)         
            fp =-lambda*( T(i+1,j) - T(i  ,j) )/dx;     % discrete Fourier law
            fm = 0;                                     % No heat flux (Neumann type BC)
        elseif(i == Nx)   
            fp = 0;                                     % No heat flux (Neumann type BC)
            fm =-lambda*( T(i  ,j) - T(i-1,j) )/dx;     % discrete Fourier law               
        else               
            fp =-lambda*( T(i+1,j) - T(i  ,j) )/dx;     % discrete Fourier law
            fm =-lambda*( T(i  ,j) - T(i-1,j) )/dx;     % discrete Fourier law
        end

        % the swipe in y direction
        if (j == 1)        
            gp =-lambda*( T(i,j+1) - T(i  ,j) )/dy;    % discrete Fourier law
            gm = 0;                                    % No heat flux (Neumann type BC) 
        elseif(j == Ny)   
            gp = 0;                                    % No heat flux (Neumann type BC)
            gm =-lambda*( T(i,j) - T(i,j-1) )/dy;      % discrete Fourier law               
        else                
            gp =-lambda*( T(i,j+1) - T(i,j  ) )/dy;    % discrete Fourier law
            gm =-lambda*( T(i,j  ) - T(i,j-1) )/dy;    % discrete Fourier law
        end

        AT(i,j) = T(i,j) + dt/dx*(fp - fm) + dt/dy*(gp - gm);  % Finite Volume update
    end
end


end