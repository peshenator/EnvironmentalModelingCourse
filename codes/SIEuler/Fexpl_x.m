% explicit (only advection terms) flux for the 1D Euler equations
% q = [rho, rho*u, rho*E]
function f = Fexpl_x(q)

    f = zeros(size(q));
    f(1,:) = q(2,:);                    % rho*u
    f(2,:) = q(2,:).^2./q(1,:);         % rho*u^2
    f(3,:) = 0.5*q(2,:).^3./q(1,:).^2;  % u*rhoK

end