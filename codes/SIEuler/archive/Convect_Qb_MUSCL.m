% Advect the cell-centerd quantities (only explicit part) using the MUSCL
% scheme with the Rusanov or FORCE numerical flux
function [Qbnew,Qxnew] = Convect_Qb_MUSCL(Qb,Qx,dt,dx)
    global Nx QL QR numFlux MUSCL;
    
    nVar = (size(Qx,1));
    dtdx = dt/dx;
    dxdt = dx/dt;
    if (MUSCL == 1)
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
    else
        Qmx = Qb;
        Qpx = Qb;
    end

    % compute PHYSical fluxes (at the barycenters x_i)
    fmb = Fexpl_x(Qmx);
    fpb = Fexpl_x(Qpx);
    
    fmx = zeros(nVar,Nx);
    fpx = zeros(nVar,Nx);
    if (numFlux == 1)
    %% ******** option1: Rusanov numerical flux at x_{i+1/2}
    umx = Qx(2,:)./Qx(1,:);
    % use ghost cells at the boundaries:
        fmx(:,1 ) = 0.5*( fmb(:,1)    + Fexpl_x(QL) - abs(umx(1   )).*( Qmx(:,1) - QL  ) );   % left  BC
        fpx(:,Nx) = 0.5*( Fexpl_x(QR) + fpb(:,Nx)   - abs(umx(Nx+1)).*( QR - Qpx(:,Nx) ) );   % right BC
        % Rusanov numerical flux at x_{i+1/2}
        for k = 1:nVar
            fmx(k,2:Nx  ) = 0.5*( fmb(k,2:Nx) + fpb(k,1:Nx-1) - abs(umx(2:Nx)).*( Qmx(k,2:Nx) - Qpx(k,1:Nx-1) ) );
            fpx(k,1:Nx-1) = fmx(k,2:Nx);
        end
    else
    %% ******** option2: FORCE numerical flux (average of Lax-Wendrof and Lax-Friedrichs):
        % Left BC:
        Q_LW = 0.5*( Qmx(:,1) + QL ) - 0.5*dtdx*( fmb(:,1) - Fexpl_x(QL) );
        F_LW = Fexpl_x(Q_LW);
        F_LF = 0.5*( fmb(:,1) + Fexpl_x(QL) ) - 0.5*dxdt*( Qmx(:,1) - QL );
        fmx(:,1 ) = 0.5*(F_LW + F_LF);
        
        % inner cells:
        Q_LW = 0.5*( Qmx(:,2:Nx) + Qpx(:,1:Nx-1) ) - 0.5*dtdx*( fmb(:,2:Nx) - fpb(:,1:Nx-1) );
        F_LW = Fexpl_x(Q_LW);
        F_LF = 0.5*( fmb(:,2:Nx) + fpb(:,1:Nx-1) ) - 0.5*dxdt*( Qmx(:,2:Nx) - Qpx(:,1:Nx-1) );
        fmx(:,2:Nx  ) = 0.5*(F_LW + F_LF);
        fpx(:,1:Nx-1) = fmx(:,2:Nx);
        
        % right BC:
        Q_LW = 0.5*( QR + Qpx(:,Nx) ) - 0.5*dtdx*( Fexpl_x(QR) - fpb(:,Nx) );
        F_LW = Fexpl_x(Q_LW);
        F_LF = 0.5*( Fexpl_x(QR) + fpb(:,Nx) ) - 0.5*dxdt*( QR - Qpx(:,Nx) );
        fpx(:,Nx) = 0.5*(F_LW + F_LF);
    end

    Qbnew = Qb - dtdx*( fpx - fmx );
    
    % avarage to the x-faces:
    Qxnew = zeros(nVar,Nx+1);
    Qxnew(:,1   ) = QL;
    Qxnew(:,2:Nx) = 0.5*(Qbnew(:,2:Nx) + Qbnew(:,1:Nx-1));
    Qxnew(:,Nx+1) = QR;


    end

                % [4] FORCE flux
                % qLW = 0.5*(q(i)+q(i+1)) - 0.5*dt/dx*(f(q(i+1))-f(q(i)));
                % fLWp = f(qLW); 
                % fLFp = 0.5*(f(q(i))+f(q(i+1))) - 0.5*dx/dt*(q(i+1)-q(i));
                % fp = 0.5*(fLWp+fLFp); % f_{i+1/2}
                % qLW = 0.5*(q(i-1)+q(i)) - 0.5*dt/dx*(f(q(i))-f(q(i-1)));
                % fLWm =  f(qLW); 
                % fLFm = 0.5*(f(q(i-1))+f(q(i))) - 0.5*dx/dt*(q(i)-q(i-1));
                % fm = 0.5 * (fLWm+fLFm); % f_{i-1/2}