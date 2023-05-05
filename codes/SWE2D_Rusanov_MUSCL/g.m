% Flux in y direction
function  flux = g(Q)
global gravity
flux(1,1) = Q(3);
flux(2,1) = Q(2)*Q(3)/Q(1);
flux(3,1) = Q(3)^2/Q(1) + 0.5*gravity*Q(1)^2;
end