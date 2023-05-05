% Numerical flux in x direction
function flux = RusanovSW_x(QL,QR,FL,FR,sL,sR)
    smax = max( max(abs(sL)),max(abs(sR)) );
    flux = 0.5*( FR + FL ) - 0.5*smax*(QR - QL) ;
end