function qfy = c2fy(qc)
    [Nx,Ny] = size(qc);
    qfy = zeros(Nx,Ny+1); 
    qfy(:,2:Ny) = 0.5*(qc(:,2:end) + qc(:,1:end-1));
    qfy(:   ,1) = 0.5*(qc(:,1) + qc(:,Ny));
    qfy(:,Nx+1) = 0.5*(qc(:,Ny) + qc(:,1));
end