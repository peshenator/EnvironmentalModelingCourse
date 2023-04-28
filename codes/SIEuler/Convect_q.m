function [qnew, qxnew] = Convect_q(q,ux,ivar)
global Nx dt dx QL QR;

qxnew = zeros(size(ux));
f = zeros(size(ux));

dtdx = dt/dx;

f(1)    = QL(2)/QL(1)*QL(ivar);
f(2:Nx) = 0.5*( ux(2:Nx).*( q(2:Nx) + q(1:Nx-1) ) - abs(ux(1:Nx-1)).*( q(2:Nx) - q(1:Nx-1) ) );
f(Nx+1) = QR(2)/QR(1)*QR(ivar);

qnew = q - dtdx*( f(2:Nx+1) - f(1:Nx) );    % q at the cell-centers

% average q to the staggered locations
qxnew(1)    = QL(ivar);
qxnew(2:Nx) = 0.5*( qnew(2:Nx) + qnew(1:Nx-1) );  
qxnew(Nx+1) = QR(ivar);

end