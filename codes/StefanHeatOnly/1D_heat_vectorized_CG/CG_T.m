% Matrix-free conjugate gradient for A*x = rhs
function [x,k,err] = CG_T(rhs,Km,Kp,dQ12)

tol = 1e-12;

N = numel(rhs);
x = rhs;                     % initial guess for the solution
r = rhs - MatVecProd_T(x,Km,Kp,dQ12);    % initial residual
p = r;                     % first vector of the basis of conjugate vectors
err = sum(sum(r.*r));           % error

for k=1:N
    if (err < tol )
        return
    end
    Ap = MatVecProd_T(p,Km,Kp,dQ12);
    alpha = err/sum(sum(p.*Ap));    % coefficients in x = alpha_i*p_i
    x = x + alpha*p;                % next approximation to the solution
    r = r - alpha*Ap;               % update the residual
    err_new = sum(sum(r.*r));       % new error
    p = r + err_new/err*p;          % the new conjugate direction
    err = err_new;
end

% disp(strcat('Conjugate Gradient does not converge, res = ',num2str(err)))

end