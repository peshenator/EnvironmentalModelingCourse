function qav = avx(q)
    % computes the average of a array q in the x-direction
    
    qav = 0.5*( q(2:end,:) + q(1:end-1,:) );

end