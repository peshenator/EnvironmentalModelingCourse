%  function for periodic BC on the top (right border for y)
function ind = bcj(j,J)
    if j==J
        ind = 1;
    else
        ind = j;
    end
end