% Derivative of the constitutive flux for the bathymetry PDE

function [fx,fy] = Diffbflux(u,v,H)
    global sA sm sUc sphi;
        % power-law DOI:10.1016/j.advwatres.2009.02.006
    
    
        fx = 1/(1-sphi)*sA*sm*(u - sUc).^(sm-1);
        fy = 1/(1-sphi)*sA*sm*(v - sUc).^(sm-1);
    end