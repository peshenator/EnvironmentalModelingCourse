function b = bathymetry(x,y)
% Bathymetry function for the 2D shallow water equations
b = zeros(size(x,1),size(y,2));

b =-0.2*(sqrt(x.^2 + y.^2) < 1);

end