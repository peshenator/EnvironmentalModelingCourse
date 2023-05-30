% The preconditioned conjugate gradient method to solve
% a linear systems of the type Ax = b
% The preconditioned systen is: E^{-1}*A*E^{-t}*y = E^{-1}*b
% where y = E^t*x and M = E*E^t
%
% Input : known vector b
% Output: x, the solution of the system
%
% The user must supply the product A*x in a function
% called "MatVecProd"
% and the precondioner matrix M^{-1}, matri M must be M = M^t > 0 and
% fixed, i.e. not changing during the iterations.
%
% Source: https://en.wikipedia.org/wiki/Conjugate_gradient_method#The_preconditioned_conjugate_gradient_method
function [x,err,k] = CGsolverPreCond(b,MatVecProd) 

    x   = b;              % initial guess
    N   = numel(b);       % 
    tol = 1e-14;          %  
    r   = b - MatVecProd(x); % initial residual 
    z   = PreCond(r); 
    p   = z;  
    for k=1:N
        err = sqrt(sum(sum(r.*r))); 
        if(err<tol)
            return
        end
        Ap     = MatVecProd(p);
        tmp    = sum(sum(r.*z));
        alpha  = tmp/sum(sum(p.*Ap));
        x      = x + alpha*p;
        r      = r - alpha*Ap;
        z      = PreCond(r);
        beta   = sum(sum(z.*r))/tmp;
        p      = z + beta*p;
    end
    
    disp(strcat('CG does NOT converge, residual = ',num2str(err)))

end