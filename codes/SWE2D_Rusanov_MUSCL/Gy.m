% Flux in y direction
function  flux = Gy(Q)
global gravity
flux = zeros(size(Q));
flux(1,:,:) = Q(3,:,:);
flux(2,:,:) = Q(2,:,:).*Q(3,:,:)./Q(1,:,:);
flux(3,:,:) = Q(3,:,:).^2./Q(1,:,:) + 0.5*gravity*Q(1,:,:).^2;
end