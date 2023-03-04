% Matrix-free conjugate gradient for A*x = b
function x = ConjGradOpt_T(b)

tol = 1e-12;

N = numel(b);
x = b;                     % initial guess for the solution
r = b - MultOpt2D_T(x);    % initial residual
p = r;                     % first vector of the basis of conjugate vectors
err = sum(sum(r.*r));           % error

for k=1:N
    if (err < tol )
        return
    end
    Ap = MultOpt2D_T(p);
    alpha = err/sum(sum(p.*Ap)); % coefficients in x = alpha_i*p_i
    x = x + alpha*p;        % next approximation to the solution
    r = r - alpha*Ap;       % update the residual
    err_new = sum(sum(r.*r));    % new error
    p = r + err_new/err*p;  % the new conjugate direction
    err = err_new;
end

disp(strcat('Conjugate Gradient does not converge, res = ',num2str(err)))

end