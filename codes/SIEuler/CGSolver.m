% Matrix-free conjugate gradient for A*x = b
function [x,k,err] = CGSolver(rhs,hx,dt,dx)

tol = 1e-14;

N = numel(rhs);
x = rhs;                    % initial guess for the solution
r = rhs - MatVecProdCG(x,hx,dt,dx);    % initial residual
p = r;                      % first vector of the basis of conjugate vectors
err = sum(sum(r.*r));       % error
err0 = err;

for k=1:N
    if (sqrt(err) < tol) % tol*sqrt(err0)
        return
    end
    Ap = MatVecProdCG(p,hx,dt,dx);
    alpha = err/sum(sum(p.*Ap));    % coefficients in x = alpha_i*p_i
    x = x + alpha*p;                % next approximation to the solution
    r = r - alpha*Ap;               % update the residual
    err_new = sum(sum(r.*r));       % new error
    p = r + err_new/err*p;          % the new conjugate direction
    err = err_new;
end

disp(strcat('CG does NOT converge, residual = ',num2str(err)))

end