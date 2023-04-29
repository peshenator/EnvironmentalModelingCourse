function rhoEx = rhoE_x(qx,px)
global gam;

% internal + kinetic
rhoEx = px/(gam-1) + 0.5*qx(2,:).^2./qx(1,:);


