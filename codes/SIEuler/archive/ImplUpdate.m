% Here, we solve the pressure subsystem with the CG method (see the paper DOI: 10.1016/j.jcp.2017.03.038)
function [qbnew,qxnew] = ImplUpdate(qb0,qx0,qb,qx,dt,dx)
global Nx QL QR pL pR gam;

qxnew = qx;

pb = pressure(qb(1,:),qb(2,:)./qb(1,:),qb(3,:));

px = zeros(1,Nx+1);
hx = zeros(1,Nx+1);
% first compute average p to x_{i+1/2}
px(1   ) = pL;
px(2:Nx) = 0.5*( pb(2:Nx) + pb(1:Nx-1) );
px(Nx+1) = pR;
%
hx(1   ) = pL/((gam-1)*QL(1)) + pL/QL(1);
hx(2:Nx) = px(2:Nx)./((gam-1)*qx(1,2:Nx)) + px(2:Nx)./qx(1,2:Nx);
hx(Nx+1) = pR/((gam-1)*QR(1)) + pR/QR(1);
%

% update m**:
qxnew(2,1   ) = QL(2);
qxnew(2,2:Nx) = qx0(2,2:Nx) - dt/dx*( pb(2:Nx) - pb(1:Nx-1) );
qxnew(2,Nx+1) = QR(2);

% update barycenter quantities:
qbnew(1:2,:) = 0.5*(qxnew(1:2,2:Nx+1) + qxnew(1:2,1:Nx));
% update rhoEb**:
qbnew(3,:) = qb0(3,:) - dt/dx*( hx(2:Nx+1).*qx(2,2:Nx+1) - hx(1:Nx).*qx(2,1:Nx) );

end