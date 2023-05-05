% Matlab function for the matrix-free conjugate gradient method 
% which solves the linear system A*x = b, 
% A is supposed to be symmetric and positive definite 
function [x,err,k] = CGsolver(b,MatVectProd,Lxm,Lxp,Lym,Lyp)
N = numel(b);       % get the size of the system 
x = b;              % initial guess 
r = b-MatVectProd(x,Lxm,Lxp,Lym,Lyp);    % initial residual = steepest descent direction 
p = r;                  % initial search direction = steepest descent direction 
err = sum(sum(r.*r));  % compute the square of r

tol = 1e-12;            % set a relative tolerance 



for k=1:N
    if(sqrt(err)<tol) 
        % we have reached the expected relative tolerance, hence return to the calling function 
        return 
    end 
    Ap = MatVectProd(p,Lxm,Lxp,Lym,Lyp);            % compute the matrix-vector product 
    alpha = err/sum(sum(p.*Ap));    % alpha minimizes the 1D problem 
    x = x + alpha*p;                % update the solution 
    r = r - alpha*Ap;               % update the residual r_k+1 = b-A*x_k+1 
    errnew = sum(sum(r.*r));        % compute the square of the new r 
    p = r + errnew/err*p;           % the new search direction (conjugate w.r.t. all previous directions)
    err = errnew;                   % overwrite the old err with the new one 
end
% if we have arrived here, we have not reached the tolerance 
disp(sprintf(' CG method did not converge. res = %e', sqrt(err) ))



















