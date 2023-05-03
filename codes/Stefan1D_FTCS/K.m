function k = K(T)

    global Tc KL KS 

    % solid
    k = KS*(T < Tc);
    
    % liquid
    k = k + KL*(T >= Tc);