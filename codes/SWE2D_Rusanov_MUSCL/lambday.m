% Eigenvalues of the Jacobian matrix of the flux in y direction
% A2 = dg/dQ
function L = lambday(Q)
global gravity
L = zeros(size(Q));
v = Q(3,:,:)./Q(1,:,:);
c = sqrt(Q(1,:,:)*gravity);
L(1,:,:) = v - c;
L(2,:,:) = v; 
L(3,:,:) = v + c;
end