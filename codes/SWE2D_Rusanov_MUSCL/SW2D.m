% Solver of shallow water equations in 2D (cartesian grids)
clear
clc
global gravity;
gravity = 9.81;
%% Computational domain
xL = -1; xR = 1;
yL = -1; yR =  1;
time = 0;
tend = 0.05;
%% Mesh generation
Nx = 50;
Ny = 50;
CFL = 0.9;
dx = (xR-xL)/Nx;
dy = (yR-yL)/Ny;
for i=1:Nx
    x(i) = xL +dx/2 +(i-1)*dx;
end
for j=1:Ny
    y(j) = yL +dy/2 +(j-1)*dy;
end
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

%% FV solver
NMAX = 1000;
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
    for i=1:Nx
        % Flux in x direction
        for j=1:Ny
            if (i==1)
                QBC = Q(:,i,j);
                QBC(2) = - QBC(2);
                fp = RusanovSW_x(Q(:,i,j),Q(:,i+1,j));
                fm = RusanovSW_x(QBC,Q(:,i,j));
            elseif (i==Nx)
                QBC = Q(:,i,j);
                QBC(2) = - QBC(2);
                fp = RusanovSW_x(Q(:,i,j),QBC);
                fm = RusanovSW_x(Q(:,i-1,j),Q(:,i,j));
            else
                fp = RusanovSW_x(Q(:,i,j),Q(:,i+1,j));
                fm = RusanovSW_x(Q(:,i-1,j),Q(:,i,j));
            end
            % Flux in y direction
            if (j==1)
                QBC = Q(:,i,j);
                QBC(3) = - QBC(3);
                gp = RusanovSW_y(Q(:,i,j),Q(:,i,j+1));
                gm = RusanovSW_y(QBC,Q(:,i,j));
            elseif (j==Ny)
                QBC = Q(:,i,j);
                QBC(3) = - QBC(3);
                gp = RusanovSW_y(Q(:,i,j),QBC);
                gm = RusanovSW_y(Q(:,i,j-1),Q(:,i,j));
            else
                gp = RusanovSW_y(Q(:,i,j),Q(:,i,j+1));
                gm = RusanovSW_y(Q(:,i,j-1),Q(:,i,j));
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
    subplot(2,1,1)
    surf(x,y,squeeze(Q(1,:,:))','EdgeColor','none','FaceColor','interp')
    axis([xL xR yL yR 0 2])
    title(strcat('Time = ',num2str(time)))

    subplot(2,1,2)
    plot(x,Q(1,:,Ny/2)','LineWidth',1.25)
    grid on
    xlim([xL xR])
    pause(0.001)
end


