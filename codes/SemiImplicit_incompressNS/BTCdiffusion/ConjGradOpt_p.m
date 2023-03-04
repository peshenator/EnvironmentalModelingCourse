% Conjugate Gradient MATRIX-FREE method for the linear system A*x = b
% INPUT:
%     b -- right hand-side 
function [x,err,i] = CGopt2D(b)

tol = 1e-12;

N = numel(b);           % get the size of the system
x = b;                  % initial guess
r = b - MultOpt2D_p(x); % initial residual 
p = r;                  % initiate p_1
err = sum(sum(r.*r));   % compute the norm of the residual = r(1)^2 + r(2)^2 +....
err0 = err;             % save the very first residual

 for i=1:N
     if (sqrt(err) < tol*sqrt(err0))
         return         % we have reached the tolerance, and we can get back to the call programm
     end
     
     Ap = MultOpt2D_p(p);          % Matrix-vectors product A*p_k
     alpha = err/sum(sum(p.*Ap));  % compute alpha_k    
     x = x + alpha*p;
     r = r - alpha*Ap;             % update the residual, it is r_{k+1} in the lecture 18
     errnew = sum(sum(r.*r));      % new error
     p = r + errnew/err*p;         % new conjugate vector
     err = errnew;
 end
 disp(strcat('Conjugate Gradient did not converge, residual = ',num2str(sqrt(err))));