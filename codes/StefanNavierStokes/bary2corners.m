function qc = bary2corners(q,qL,qR,qB,qT)

% bary2corners - convert barycentric values to corner values
% qL, qR, qB, qT are the values at the left, right, bottom, and top
% corners of the cell, respectively.  q is the barycentric value.

[Nx,Ny] = size(q);
qc = zeros(Nx+1,Ny+1);

qc(2:Nx,1:Ny) = 0.5*( q(2:Nx,:   ) + q(1:Nx-1,:    ));
qc(2:Nx,2:Ny) = 0.5*(qc(2:Nx,2:Ny) + qc(2:Nx,1:Ny-1));
qc(1     ,1:Ny) = qL;
qc(Nx+1  ,1:Ny) = qR;
qc(1:Nx  ,1   ) = qB;
qc(1:Nx+1,Ny+1) = qT;

end