% Advect the cell-centerd quantities (only explicit part)
function [Qbnew,Qxnew] = Convect_Qb_MUSCL(Qb,Qx,dt,dx)
    global Nx QL QR;
    
    nVar = (size(Qx,1));
    dtdx = dt/dx;

    %% Piecewise linear reconstruction in space and time
    slopex = zeros(nVar,Nx);
    slopex(:,2:Nx-1) = minmodarray(Qb(:,2:Nx-1)-Qb(:,1:Nx-2),Qb(:,3:Nx)-Qb(:,2:Nx-1));
  
    % extrapolating Q to the x-faces
    Qmx = Qb - 0.5*slopex; % extrapolated data at (i-1/2,j)
    Qpx = Qb + 0.5*slopex; % extrapolated data at (i+1/2,j)
    
    % if doing an IMEX, the extrapolation in time should be OFF
    Qt =-( Fexpl_x(Qpx) - Fexpl_x(Qmx) )/dx;
    % update Qmx=Q_{i-1/2,j} and Qpx=Q_{i+1/2,j} in time by dt/2
    Qmx = Qmx + 0.5*dt*Qt;
    Qpx = Qpx + 0.5*dt*Qt;

    % compute PHYSical fluxes (at the barycenters x_i)
    fmb = Fexpl_x(Qmx);
    fpb = Fexpl_x(Qpx);
    
    % compute NUMerical fluxes (at the interfaces x_{i+1/2})
    umx = Qmx(2,:)./Qmx(1,:);
    upx = Qpx(2,:)./Qpx(1,:);

    fmx = zeros(nVar,Nx);
    fpx = zeros(nVar,Nx);
    fmx(:,1 ) = Fexpl_x(QL);
    fpx(:,Nx) = Fexpl_x(QR);
    for k = 1:nVar
        fmx(k,2:Nx  ) = 0.5*( fmb(k,2:Nx) + fpb(k,1:Nx-1) - abs(umx(2:Nx)).*( Qmx(k,2:Nx) - Qpx(k,1:Nx-1) ) );
        fpx(k,1:Nx-1) = fmx(k,2:Nx);
    end

    Qbnew = Qb - dtdx*( fpx - fmx );
    
    % avarage to the x-faces:
    Qxnew = zeros(nVar,Nx+1);
    Qxnew(:,1   ) = QL;
    Qxnew(:,2:Nx) = 0.5*(Qbnew(:,2:Nx) + Qbnew(:,1:Nx-1));
    Qxnew(:,Nx+1) = QR;


    end