% Initial data for the Riemann problems for the ideal gas equation of state
function [rhoL, uL, pL, rhoR, uR, pR, xL, xR, tend] = RPinit(nRP)

if (nRP == 0 )
    rhoL = 1e3;     % density left
    uL   = 1;       % u left
    pL   = 1e5;     % pressure left

    rhoR = 0.01;    % density right
    uR   = 1;       % u right
    pR   = 1e5;     % pressure right
    
    xL =-0.5;
    xR = 0.5;
    tend = 0.2;
elseif(nRP == 1 )
    rhoL = 1;       % density left
    uL   = 0;       % u left
    pL   = 1;       % pressure left

    rhoR = 0.125;   % density right
    pR   = 0.1;     % pressure right
    uR   = 0;       % u right
    
    xL =-0.5;
    xR = 0.5;
    tend = 0.2;
elseif( nRP == 2 )
    rhoL = 0.445;   % density left
    uL   = 0.698;   % u left
    pL   = 3.528;   % pressure left

    rhoR = 0.5;     % density right
    uR   = 0;       % u right
    pR   = 0.571;   % pressure right

    xL =-0.5;
    xR = 0.5;
    tend = 0.14;
elseif( nRP == 3 )
    rhoL = 1;       % density left
    uL   = 0;       % u left
    pL   = 1e3;     % pressure left

    rhoR = 1;       % density right
    uR   = 0;       % u right
    pR   = 0.1;     % pressure right
    
    xL =-0.5;
    xR = 0.5;
    tend = 0.012;
elseif( nRP == 4 )
    rhoL = 5.99924; % density left
    uL   = 19.5975; % u left
    pL   = 460.894; % pressure left

    rhoR = 5.99924;  % density right
    uR   =-6.19633;  % u right
    pR   = 46.095;   % pressure right
    
    xL =-1;
    xR = 1;
    tend = 0.035;
elseif( nRP == 5 )
    rhoL = 1;       % density left
    uL   =-1;       % u left
    pL   = 0.4;     % pressure left

    rhoR = 1;       % density right
    uR   = 1;       % u right
    pR   = 0.4;     % pressure right
    
    xL =-0.5;
    xR = 0.5;
    tend = 0.15;
elseif( nRP == 6 )
    rhoL = 1;       % density left
    uL   = 2;       % u left
    pL   = 0.1;     % pressure left

    rhoR = 1;       % density right
    uR   =-2;       % u right
    pR   = 0.1;     % pressure right
    
    xL =-0.5;
    xR = 0.5;
    tend = 0.8;
elseif( nRP == 7 )
    rhoL = 1;       % density left
    uL   = 0.75;    % u left
    pL   = 1;       % pressure left

    rhoR = 0.125;   % density right
    uR   = 0;       % u right
    pR   = 0.1;     % pressure right
    
    xL =-0.5;
    xR = 0.5;
    tend = 0.2;
elseif( nRP == 8 )
    rhoL = 1;       % density left
    uL   = 0.1;     % u left
    pL   = 0.1;       % pressure left

    rhoR = 1;       % density right
    uR   =-0.1;     % u right
    pR   = 0.1;       % pressure right
    
    xL =-0.5;
    xR = 0.5;
    tend = 0.4;
end
end
