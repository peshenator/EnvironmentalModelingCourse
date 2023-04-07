function qav = avy(q)
    % compute the average of a array q in the y-direction
    
    qav = 0.5*( q(:,2:end) + q(:,1:end-1) );

end