% Coompute M^{-1}*b = M^{-1}*Ax, where M = M^t > 0 is the preconditioner
% for the linear system A*x = b.
% P = M^{-1}
function PAx = PreCond(Ax) 
global a b c di Nx Nz 

PAx = zeros(Nz+1,Nx);
for i=1:Nx
    PAx(:,i) = Thomas(a(:,i),b(:,i) + di(:,i),c(:,i),Ax(:,i)); 
end
