% Flux in x direction
function  flux = Fx(Q)
global gravity
flux = zeros(size(Q));
flux(1,:,:) = Q(2,:,:);
flux(2,:,:) = Q(2,:,:).^2./Q(1,:,:) + 0.5*gravity*Q(1,:,:).^2;
flux(3,:,:) = Q(2,:,:).*Q(3,:,:)./Q(1,:,:);
end