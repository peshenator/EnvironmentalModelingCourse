% The conjugate gradient method to solve
% linear systems of the type Ax = b
% Input: known vector b
% Output: x, the solution of the system
% The user must supply the product A*x in a function
% called "MatVecProd"
function [x,res,k] = CGsolverPreCond(b,MatVecProd) 
x = b;              % initial guess
N = numel(b);       % 
tol = 1e-14;        %  
r = b - MatVecProd(x); % initial residual 
z = PreCond(r); 
p = z;  
for k=1:2*N
    res = sqrt(sum(sum(sum(r.*r)))); 
    if(res<tol)
        % abbiamo gia trovato la soluzione 
        return
    end
    v = MatVecProd(p);
    alpha = sum(sum(sum(r.*z)));
    lambda = alpha/sum(sum(sum(p.*v)));
    x = x + lambda*p;
    r = r - lambda*v;
    z = PreCond(r);
    beta = sum(sum(sum(z.*r)))/alpha;
    p = z + beta*p;           % nuova dir. di ricerca 
end