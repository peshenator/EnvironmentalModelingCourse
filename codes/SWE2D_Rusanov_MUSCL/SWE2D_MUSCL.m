% 2D Finite Volume MUSCL scheme for the Shallow Water Equations
% based on the Rusanov fluxclear;
clear;
global gravity;
gravity = 9.81;
%% Computational domain
xL = -1; xR = 1;
yL = -1; yR =  1;
time = 0;
tend = 1;
%% Mesh generation
Nx = 128;
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
        r = sqrt((x(i) + 0.75)^2+(y(j) + 0.75)^2);
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
% Time loop
for n=1:10000000000
    % Computation of the time step 
    sx = lambdax(Q);
    sy = lambday(Q);
    amax_x = max( max( max( abs( sx ) ) ) );
    amax_y = max( max( max( abs( sy ) ) ) );
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
    % X - direction
    slopeX(:,2:Nx-1,:) = minmodarray(Q(:,2:Nx-1,:)-Q(:,1:Nx-2,:),Q(:,3:Nx,:)-Q(:,2:Nx-1,:));
    %
    % Y - direction
    slopeY(:,:,2:Ny-1) = minmodarray(Q(:,:,2:Ny-1)-Q(:,:,1:Ny-2),Q(:,:,3:Ny)-Q(:,:,2:Ny-1));
    % Compute extrapolated X-data at the cell boundaries
    Qxm = Q - 0.5*slopeX;
    Qxp = Q + 0.5*slopeX;
    % Compute extrapolated Y-data at the cell boundaries
    Qym = Q - 0.5*slopeY;
    Qyp = Q + 0.5*slopeY;
    % Compute time derivative from the PDE:
    Q_t =-(Fx(Qxp) - Fx(Qxm))/dx - (Gy(Qyp) - Gy(Qym))/dy;
    % Evolve in time boundary extrapolated X-data
    Qxm = Qxm + 0.5*dt*Q_t;
    Qxp = Qxp + 0.5*dt*Q_t;
    % Evolve in time boundary extrapolated Y-data
    Qym = Qym + 0.5*dt*Q_t;
    Qyp = Qyp + 0.5*dt*Q_t;
 
    % update characteristic speeds at the xtrapolated data
    sxm = lambdax(Qxm); sxp = lambdax(Qxp);
    sym = lambday(Qym); syp = lambday(Qyp);
    % update the fluxes at the extrapolated data
    fxm = Fx(Qxm);  fxp = Fx(Qxp);
    gym = Gy(Qym);  gyp = Gy(Qyp);

%% Standard first-order finite volume
    for i=1:Nx
        % Flux in x direction
        for j=1:Ny
            if (i==1)
                QBC    = Q(:,i,j);
                QBC(2) =-QBC(2);
                fm = Rusanov(QBC       ,Qxm(:,i  ,j),Fx(QBC)     ,fxm(:,i  ,j),lambdax(QBC),sxm(:,i  ,j));
                fp = Rusanov(Qxp(:,i,j),Qxm(:,i+1,j),fxp(:,i,j  ),fxm(:,i+1,j),sxp(:,i,j)  ,sxm(:,i+1,j));
            elseif (i==Nx)
                QBC    = Q(:,i,j);
                QBC(2) =-QBC(2);
                fm = Rusanov(Qxp(:,i-1,j),Qxm(:,i,j),fxp(:,i-1,j),fxm(:,i  ,j),sxp(:,i-1,j),sxm(:,i  ,j));
                fp = Rusanov(Qxp(:,i,j  ),QBC       ,fxp(:,i,j  ),Fx(QBC)     ,sxp(:,i  ,j),lambdax(QBC));
            else
                fm = Rusanov(Qxp(:,i-1,j),Qxm(:,i  ,j),fxp(:,i-1,j),fxm(:,i  ,j),sxp(:,i-1,j),sxm(:,i  ,j));
                fp = Rusanov(Qxp(:,i,j  ),Qxm(:,i+1,j),fxp(:,i,j  ),fxm(:,i+1,j),sxp(:,i  ,j),sxm(:,i+1,j));          
            end
            % Flux in y direction
            if (j==1)
                QBC    = Q(:,i,j);
                QBC(3) =-QBC(3);
                gm = Rusanov(QBC       ,Qym(:,i,j)  ,Gy(QBC)     ,gym(:,i,j  ),lambday(QBC),sym(:,i,j  ));
                gp = Rusanov(Qyp(:,i,j),Qym(:,i,j+1),gyp(:,i,j  ),gym(:,i,j+1),syp(:,i,j  ),sym(:,i,j+1));
            elseif (j==Ny)
                QBC    = Q(:,i,j);
                QBC(3) =-QBC(3);
                gm = Rusanov(Qyp(:,i,j-1),Qym(:,i,j),gyp(:,i,j-1),gym(:,i,j  ),syp(:,i,j-1),sym(:,i,j  ));
                gp = Rusanov(Qyp(:,i,j),QBC         ,gyp(:,i,j  ),Gy(QBC)     ,syp(:,i,j  ),lambday(QBC));
            else
                gm = Rusanov(Qyp(:,i,j-1),Qym(:,i,j  ),gyp(:,i,j-1),gym(:,i,j  ),syp(:,i,j-1),sym(:,i,j  ));
                gp = Rusanov(Qyp(:,i,j  ),Qym(:,i,j+1),gyp(:,i,j  ),gym(:,i,j+1),syp(:,i,j  ),sym(:,i,j+1));
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


