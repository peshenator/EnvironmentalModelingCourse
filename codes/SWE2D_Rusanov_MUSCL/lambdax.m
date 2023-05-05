% Eigenvalues of the Jacobian matrix of the flux in x direction
% A1 = df/dQ
function L = lambdax(Q)
global gravity
L = zeros(size(Q));
u = Q(2,:,:)./Q(1,:,:);
c = sqrt(Q(1,:,:)*gravity);
L(1,:,:) = u - c;
L(2,:,:) = u;
L(3,:,:) = u + c;
end