function gamma=Newton(gamma0)
gamma = gamma0;      % initial guess 
tol   = 1e-11;       % tolerance 
MaxNewton = 100;     % maximum number of Newton iterations 
for iNewton=1:MaxNewton
    gk = g(gamma);   % evaluate the function 
    res = abs(gk);   % compute the residual (error) 
    disp(sprintf('iNewton = %d, res = %e',iNewton,res)) 
    if(res<tol)
        break
    end
    dgk = dg(gamma); % compute the derivative of g 
    dgamma = -gk/dgk;  
    gamma = gamma + dgamma; 
end