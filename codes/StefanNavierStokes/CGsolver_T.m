% Matrix-free conjugate gradient for A*x = b
function [x,k,err] = CGsolver_T(rhs,MatVecProd,Kmx,Kpx,Kmy,Kpy,dQ12)

tol = 1e-12;

N = numel(rhs);
x = rhs;                     % initial guess for the solution
r = rhs - MatVecProd(x,Kmx,Kpx,Kmy,Kpy,dQ12);    % initial residual
p = r;                     % first vector of the basis of conjugate vectors
err = sum(sum(r.*r));           % error

for k=1:N
    if (err < tol )
        return
    end
    Ap = MatVecProd(p,Kmx,Kpx,Kmy,Kpy,dQ12);
    alpha = err/sum(sum(p.*Ap));    % coefficients in x = alpha_i*p_i
    x = x + alpha*p;                % next approximation to the solution
    r = r - alpha*Ap;               % update the residual
    err_new = sum(sum(r.*r));       % new error
    p = r + err_new/err*p;          % the new conjugate direction
    err = err_new;
end

disp(strcat('CG T does not converge, res = ',num2str(err)))

end