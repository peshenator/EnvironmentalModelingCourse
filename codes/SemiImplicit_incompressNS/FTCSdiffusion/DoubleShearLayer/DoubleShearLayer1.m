function [u,v,T] = DoubleShearLayer1(u,v,T,direct_flag)
global Nx Ny xb yb;

u1 = 0.5;
u2 =-0.5;
um = (u1-u2)/2;
L = 0.025;
d = 1e-2;
if (direct_flag == 1)
    % u
    for j = 1:Ny
        if    (yb(j) < 0.25)
            u(:,j) = u1 - um*exp(( yb(j)-0.25 )/L);
        elseif(yb(j) < 0.5)
            u(:,j) = u2 + um*exp((-yb(j)+0.25 )/L);
        elseif(yb(j)<0.75)
            u(:,j) = u2 + um*exp(( yb(j)-0.75 )/L);
        elseif(yb(j)<1.0)
            u(:,j) = u1 - um*exp((-yb(j)+0.75 )/L);
        end
    end
    % v
    for j = 1:Ny+1
        v(:,j) = d*sin(4*pi*xb);
    end
    T = u(1:Nx,:);  % passive scalar
else
    % v
    for i = 1:Nx
        if    (xb(i) < 0.25)
            v(i,:) = u1 - um*exp(( xb(i)-0.25 )/L);
        elseif(xb(i)<0.5)
            v(i,:) = u2 + um*exp((-xb(i)+0.25 )/L);
        elseif(xb(i)<0.75)
            v(i,:) = u2 + um*exp(( xb(i)-0.75 )/L);
        elseif(xb(i)<1.0)
            v(i,:) = u1 - um*exp((-xb(i)+0.75 )/L);
        end
    end
    % % u
    for i = 1:Nx+1
        u(i,:) = d*sin(4*pi*yb);
    end
    T = v(:,1:Ny);  % passive scalar
end

end