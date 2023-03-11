function qfx = c2fx(qc)
    [Nx,Ny] = size(qc);
    qfx = zeros(Nx+1,Ny); 
    qfx(2:Nx,:) = 0.5*(qc(2:end,:) + qc(1:end-1,:));
    qfx(1   ,:) = 0.5*(qc(1,:) + qc(Nx,:));
    qfx(Nx+1,:) = 0.5*(qc(Nx,:) + qc(1,:));
end