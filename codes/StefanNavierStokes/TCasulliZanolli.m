% BTCS scheme for the Stefan problem based on the nested Newton method of
% Casulli and Zanolli

function T = TCasulliZanolli(T)
    global rhoS hs Ts epsilonT
    % physical parameters of water and ice [in SI units]

    % dt = d/( max(kappaS,kappaL)/dx^2 + max(kappaS,kappaL)/dy^2 );         % time step restriction of FTCS
    % compute the coefficients of the linear part of the system on T
    [Kmx,Kpx,Kmy,Kpy,rhs] = LinearPartCoeff(T);
    T0 = T;
    
    %%  Newsted Newton method of Casulli & Zanolli
    tol = 1e-5*rhoS*hs;
    T0 = min(T0,Ts-epsilonT); % Initial guess for the outer iterations, see the paper by Casulli & Zanolli
    MaxNewton = 100;
    % ----- OUTER iterations -------
    for iouter = 1:MaxNewton
        ResOut = ResidualOuter(T0,Kmx,Kpx,Kmy,Kpy,rhs);    % Compute the residual of the outer iterations
        ResOut_norm = norm(ResOut);             % the norm of the residual
        if (ResOut_norm < tol)
            break
        end
        Talpha = T0; % we store the value of the temperature so that from now on the meaning of T0 is T^(alpha,beta)
        % ----- INNER iterations ------
        T0 = max(T0,Ts-epsilonT);   % Initial guess for the inner iterations, see the paper by Casulli & Zanolli
        for iinner = 1:MaxNewton
            [ResInner,dQ12] = ResidualInner(T0,Talpha,Kmx,Kpx,Kmy,Kpy,rhs);    % Compute the residual of the inner iterations and the Jacobian of Q=Q1-Q2
            ResInner_norm = norm(ResInner);
%             disp(strcat('        Inner iteration = ',num2str(iinner),'|| residual =',num2str(ResInner_norm)))
            if (ResInner_norm < tol)
                break
            end
            % dT = Thomas(a,b+dQ12,c,ResInner);
            [dT,CGiter,CGerr] = CGsolver_T(ResInner,@MatVecProd_T,Kmx,Kpx,Kmy,Kpy,dQ12);
            T0 = T0 - dT;
        end
    end
    T = T0;

end