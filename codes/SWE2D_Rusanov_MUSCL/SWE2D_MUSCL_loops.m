% 2D Finite Volume MUSCL scheme for the Shallow Water Equations
% based on the Rusanov fluxclear;
clear;
global gravity;
gravity = 9.81;
%% Computational domain
xL = -1; xR = 1;
yL = -1; yR =  1;
time = 0;
tend = 0.1;
%% Mesh generation
Nx = 64;
Ny = Nx;
CFL = 0.9;
dx = (xR-xL)/Nx;
dy = (yR-yL)/Ny;
x = linspace(xL + dx/2,xR - dx/2,Nx);
y = linspace(yL + dy/2,yR - dy/2,Ny);
%% Initial conditions
Q = zeros(3,Nx,Ny);
Qnew = zeros(3,Nx,Ny);
for i=1:Nx
    for j=1:Ny
        r = sqrt(x(i)^2+y(j)^2);
        if (r < 0.5)
            Q(:,i,j) = [2 ; 0 ; 0 ];
        else
            Q(:,i,j) = [0.5 ; 0 ; 0];
        end
    end
end
surf(x,y,squeeze(Q(1,:,:))')
axis([xL xR yL yR 0 2])
%% Allocate auxiliary arrays
Qxm = zeros(3,Nx,Ny);
Qxp = zeros(3,Nx,Ny);
Qym = zeros(3,Nx,Ny);
Qyp = zeros(3,Nx,Ny);
slopeX = zeros(3,Nx,Ny);
slopeY = zeros(3,Nx,Ny);
%% FV solver
NMAX = 1000000;
% Time loop
for n=1:NMAX
    % Computation of the time step 
    amax_x = 0; amax_y = 0;
    for i=1:Nx
        for j=1:Ny
            amax_x = max( amax_x, max( abs( lambdax(Q(:,i,j) ) ) ) );
            amax_y = max( amax_y, max( abs( lambday(Q(:,i,j) ) ) ) );
        end
    end
    dt = CFL / ( amax_x/dx + amax_y/dy );
    % Ensure tend
    if (time+dt>tend)
        dt = tend-time;
    end
    % Stop criterium
    if (time>=tend)
        break
    end
    
%% Piecewise linear reconstruction in space and time
    for i=1:Nx
        for j=1:Ny      
            % X - direction
            if (i==1)
                slopeX(:,i,j) = zeros(3,1);
            elseif(i==Nx)
                slopeX(:,i,j) = zeros(3,1);
            else
                slopeX(:,i,j) = minmod(Q(:,i,j)-Q(:,i-1,j),Q(:,i+1,j)-Q(:,i,j));
            end
            %
            % Y - direction
            if (j==1)
                slopeY(:,i,j) = zeros(3,1);
            elseif(j==Ny)
                slopeY(:,i,j) = zeros(3,1);
            else
                slopeY(:,i,j) = minmod(Q(:,i,j)-Q(:,i,j-1),Q(:,i,j+1)-Q(:,i,j));
            end
            % Compute extrapolated X-data at the cell boundaries
            Qxm(:,i,j) = Q(:,i,j) - 0.5*slopeX(:,i,j);
            Qxp(:,i,j) = Q(:,i,j) + 0.5*slopeX(:,i,j);
            % Compute extrapolated Y-data at the cell boundaries
            Qym(:,i,j) = Q(:,i,j) - 0.5*slopeY(:,i,j);
            Qyp(:,i,j) = Q(:,i,j) + 0.5*slopeY(:,i,j);
            % Compute time derivative from the PDE:
            Q_t =-(f(Qxp(:,i,j)) - f(Qxm(:,i,j)))/dx - ...
                (g(Qyp(:,i,j)) - g(Qym(:,i,j)))/dy;
            % Evolve in time boundary extrapolated X-data
            Qxm(:,i,j) = Qxm(:,i,j) + 0.5*dt*Q_t;
            Qxp(:,i,j) = Qxp(:,i,j) + 0.5*dt*Q_t;
            % Evolve in time boundary extrapolated Y-data
            Qym(:,i,j) = Qym(:,i,j) + 0.5*dt*Q_t;
            Qyp(:,i,j) = Qyp(:,i,j) + 0.5*dt*Q_t;
        end
    end
%% Standard first-order finite volume
    for i=1:Nx
        % Flux in x direction
        for j=1:Ny
            if (i==1)
                QBC = Q(:,i,j);
                QBC(2) = - QBC(2);
                fm = RusanovSW_x(QBC,Qxm(:,i,j));
                fp = RusanovSW_x(Qxp(:,i,j),Qxm(:,i+1,j));
            elseif (i==Nx)
                QBC = Q(:,i,j);
                QBC(2) = - QBC(2);
                fm = RusanovSW_x(Qxp(:,i-1,j),Qxm(:,i,j));
                fp = RusanovSW_x(Qxp(:,i,j),  QBC);
            else
                fm = RusanovSW_x(Qxp(:,i-1,j),Qxm(:,i,j));
                fp = RusanovSW_x(Qxp(:,i,j),  Qxm(:,i+1,j));
                if (max(abs(Q(:,i,j)-Q(:,i-1,j))) > 0)
                    Q_t = 0;
                end                 
            end
            % Flux in y direction
            if (j==1)
                QBC = Q(:,i,j);
                QBC(3) = - QBC(3);
                gm = RusanovSW_y(QBC,Qym(:,i,j));
                gp = RusanovSW_y(Qyp(:,i,j),Qym(:,i,j+1));
            elseif (j==Ny)
                QBC = Q(:,i,j);
                QBC(3) = - QBC(3);
                gm = RusanovSW_y(Qyp(:,i,j-1),Qym(:,i,j));
                gp = RusanovSW_y(Qyp(:,i,j),QBC);
            else
                gm = RusanovSW_y(Qyp(:,i,j-1),Qym(:,i,j));
                gp = RusanovSW_y(Qyp(:,i,j),Qym(:,i,j+1));
                if (max(abs(Q(:,i,j)-Q(:,i,j-1))) > 0)
                    Q_t = 0;
                end                 
            end
            % Computation of the solution at the new time step using
            % FV
            Qnew(:,i,j) = Q(:,i,j) - dt/dx*(fp-fm) - dt/dy*(gp-gm);
        end
    end        
    % Overwrite solution
    Q = Qnew;
    % Update time
    time = time +dt; 

    % Plot solution
    subplot(1,2,1)
    surf(x,y,squeeze(Q(1,:,:))','EdgeColor','none','FaceColor','interp')
    axis([xL xR yL yR 0 2])
    title(strcat('Time = ',num2str(time)))
    light('Position',[-3 -1 1],'Style','local')
    %
    subplot(1,2,2)
    plot(x,Q(1,:,Ny/2)','LineWidth',1.25)
    grid on
    xlim([xL xR])
    pause(0.001)
end


