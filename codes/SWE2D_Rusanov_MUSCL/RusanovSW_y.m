% Numerical flux in y direction
% sL, sR are the characteristic velocities at the left and right states
function flux = RusanovSW_y(QL,QR,Fl,FR,sL,sR)
smax = max( max(abs(sL)),max(abs(sR)) );
flux = 0.5*( FR + FL ) - 0.5*smax*(QR - QL) ;
end