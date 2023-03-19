function [u,v,T] = DoubleShearLayer2(u,v,T,direct_flag)
global Nx Ny yb xb;

L = 30;
d = 5e-2;

if (direct_flag == 1)
    % u
    for j = 1:Ny
        if (yb(j) <= 0.5)
            u(:,j) = tanh(L*(yb(j) - 0.25));
        else
            u(:,j) = tanh(L*(0.75 - yb(j)));
        end
    end
    % v
    for j = 1:Ny+1
        v(:,j) = d*sin(2*pi*xb);
    end
    T = u(1:Nx,:);  % passive scalar
else
    % v
    for i = 1:Nx
        if (xb(i) <= 0.5)
            v(i,:) = tanh(L*(xb(i) - 0.25));
        else
            v(i,:) = tanh(L*(0.75 - xb(i)));
        end
    end
    % % u
    for i = 1:Nx+1
        u(i,:) = d*sin(2*pi*yb);
    end
    T = v(:,1:Ny);  % passive scalar
end

end