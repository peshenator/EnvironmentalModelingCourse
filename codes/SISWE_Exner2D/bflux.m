% constitutive flux for the bathymetry PDE

function [fx,fy] = bflux(u,v,H)
global sA sm sUc sphi;
    % power-law DOI:10.1016/j.advwatres.2009.02.006

    fx = 1/(1-sphi)*sA*(u - sUc).^sm;
    fy = 1/(1-sphi)*sA*(v - sUc).^sm;
end