%  function for periodic BC
function ind = bcj_1(j,J)
    if j==0
        ind = J;
    else
        ind = j;
    end
end