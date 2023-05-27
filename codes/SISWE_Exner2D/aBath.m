% Derivative (characteristic velocity) of the constitutive flux for the bathymetry PDE

function [ax,ay] = aBath(u,v,H)
    global sA sm sUc sphi;
        % power-law DOI:10.1016/j.advwatres.2009.02.006
    
    
        ax = 1/(1-sphi)*sA*sm*(u - sUc).^(sm-1);
        ay = 1/(1-sphi)*sA*sm*(v - sUc).^(sm-1);
    end