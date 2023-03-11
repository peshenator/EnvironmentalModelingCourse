%  function for periodic BC
function ind = bcind(j,J)
    if j==0
        ind = J;
    else
        ind = j;
    end
end