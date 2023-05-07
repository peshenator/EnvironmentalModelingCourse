% Here, we solve the pressure subsystem with the CG method (see the paper DOI: 10.1016/j.jcp.2017.03.038)
function [qbnew,qxnew,res,iPicard] = PressureSubSystem(qb,qx,dt,dx)
global Nx pL pR QL QR gam MaxPicardP tolPicard;

res = 1;

hx = zeros(1,Nx+1); % enthalpy at the x-interfaces
px = zeros(1,Nx+1); % pressure at the x-interfaces

qbnew = qb;
qxnew = qx;

rhox = qx(1,:);
mx   = qx(2,:);
rhoEb= qb(3,:);

% initial guess for the pessure:
pb = pressure(qb(1,:),qb(2,:)./qb(1,:),qb(3,:));

% speed for extra dissipation
dtdx  = dt/dx;
dtdx2 = dtdx^2;

for iPicard=1:MaxPicardP
    if (res < tolPicard)
        break
    end
    %% Update the enthalpy at the interfaces x_{i+1/2}:
    % first compute average p to x_{i+1/2}
    px(1   ) = pL;
    px(2:Nx) = max(pb(2:Nx),pb(1:Nx-1));
    px(Nx+1) = pR;
    %
    hx(1   ) = pL/((gam-1)*QL(1))             + pL/QL(1);
    hx(2:Nx) = px(2:Nx)./((gam-1)*rhox(2:Nx)) + px(2:Nx)./rhox(2:Nx);
    hx(Nx+1) = pR/((gam-1)*QR(1))             + pR/QR(1);
    %
    
    %% rhs:
    % kinetic energy:
%     ux = mx./rhox;
%     rhoK = 0.25*( rhox(2:Nx+1).*ux(2:Nx+1).^2 + rhox(1:Nx).*ux(1:Nx).^2 );
    rhoK = 0.25*( mx(2:Nx+1).^2./rhox(2:Nx+1) + mx(1:Nx).^2./rhox(1:Nx) );

    rhs    = rhoEb - rhoK - dtdx*( hx(2:Nx+1).*qx(2,2:Nx+1) - hx(1:Nx).*qx(2,1:Nx) );
    rhs(1 ) = rhs(1 ) + dtdx2*hx(1 )*pL;
    rhs(Nx) = rhs(Nx) + dtdx2*hx(Nx)*pR;
    
    % solve the linear system on the pressure with the CG method:
    res = norm(MatVecProd_p(pb,hx,dt,dx) - rhs)/norm(pb);   % residual of the mildly nonlinear system
    if (res < tolPicard)
        break
    end    
    [pb,CGerr,CGk] = CGSolver(rhs,@MatVecProd_p,hx,dt,dx);
    % update m**:
    mx(1   ) = QL(2);
    mx(2:Nx) = qx(2,2:Nx) - dtdx*( pb(2:Nx) - pb(1:Nx-1) );
    mx(Nx+1) = QR(2);

end
% update face quantities (we don't care about qx(3,:)):
qxnew(2,:) = mx;

% update barycenter quantities:
qbnew(1:2,:) = 0.5*(qxnew(1:2,2:Nx+1) + qxnew(1:2,1:Nx));
% update rhoEb**:
qbnew(3,:) = rhoEb - dtdx*( hx(2:Nx+1).*mx(2:Nx+1) - hx(1:Nx).*mx(1:Nx) );


if (res > tolPicard)
    disp(strcat('Pressure system does NOT converge, res = ',num2str(res)))
elseif (isnan(res))
    disp('NaN in Pressure step')
end
end