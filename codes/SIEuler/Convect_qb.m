% Advect the cell-centerd quantities (only explicit part)
function [qbnew,qxnew] = Convect_qb(qb,qx,dt,dx)
    global Nx QL QR;
    
    nVar = (size(qx,1));
    dtdx = dt/dx;
    
    % compute PHYSical fluxes (at the barycenters x_i)
    fb = Fexpl_x(qb);
    
    % compute NUMerical fluxes (at the interfaces x_{i+1/2})
    ux = qx(2,:)./qx(1,:);
    fx = zeros(nVar,Nx+1);

    fx(:,1   ) = Fexpl_x(QL);
    fx(:,Nx+1) = Fexpl_x(QR);
    for k = 1:nVar
        fx(k,2:Nx) = 0.5*( fb(k,2:Nx) + fb(k,1:Nx-1) - abs(ux(2:Nx)).*( qb(k,2:Nx) - qb(k,1:Nx-1) ) );
    end

    qbnew = qb - dtdx*( fx(:,2:Nx+1) - fx(:,1:Nx) );

    % avarage to the x-faces:
    qxnew = zeros(nVar,Nx+1);
    qxnew(:,1   ) = QL;
    qxnew(:,2:Nx) = 0.5*(qbnew(:,2:Nx) + qbnew(:,1:Nx-1));
    qxnew(:,Nx+1) = QR;

    end