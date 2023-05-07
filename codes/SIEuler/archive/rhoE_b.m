function rhoEb = rhoE_b(qb,pb)
global gam;

% internal + kinetic
rhoEb = pb/(gam-1) + 0.5*qb(2,:).^2./qb(1,:);


