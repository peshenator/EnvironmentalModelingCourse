function [u,v,T] = BubbleAdvection(u,v,T,a)
global Nx Ny xb yb;

if norm(a) > 0
    a = a/norm(a);
end

u(:,:) = a(1);
v(:,:) = a(2);

for i=1:Nx
    for j=1:Ny
        T(i,j) = exp(-( (xb(i)-0.5)^2 + (yb(j)-0.5)^2 )/0.1^2 )^2; 
    end
end

end
