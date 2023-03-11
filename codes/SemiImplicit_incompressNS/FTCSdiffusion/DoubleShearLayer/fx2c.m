function qc = fx2c(q)
    qc = 0.5*(q(2:end,:) + q(1:end-1,:));
end