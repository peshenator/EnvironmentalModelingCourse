% matrix-free conjugate gradient method for 
% the solution of A*x = b, 
% with A symmetric positive definite
% Input : 
%   b = right hand side 
% Output: 
%   x = solution vector 
% The user has to provide a function "matop2D(x)" which 
% computes the matrix-vector product A*x 
function x=CGop2D(b)
x = b;        % initial guess 
N = numel(b); % get the problem size 
r = b-matop2D(x);    % residual 
p = r;        % initial search direction 
alpha = sum(sum(r.*r)); % square of r 
for k=1:4*N
    if(sqrt(alpha)<1e-14) 
        % norm of the residual is sufficiently small 
        break
    end
    v = matop2D(p); 
    lambda = alpha/sum(sum(p.*v));
    x = x + lambda*p; 
    r = r - lambda*v;
    alphak = alpha; 
    alpha = sum(sum(r.*r)); 
    p = r + alpha/alphak*p; 
end



