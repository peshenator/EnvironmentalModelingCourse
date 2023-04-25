function PAx = precop(Ax) 
global a b c di Nx Nz 
%PAx = Ax;
PAx = zeros(Nz+1,Nx);
for i=1:Nx
    PAx(:,i) = Thomas(a(:,i),b(:,i)+di(:,i),c(:,i),Ax(:,i)); 
end
