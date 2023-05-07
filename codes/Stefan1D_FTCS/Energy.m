function Q = Energy(T)

global Tc hL rhoS rhoL CS CL;

if (T < Tc)
    Q = rhoS*CS*(T-Tc);
else
    Q = rhoL*CL*(T-Tc) + rhoL*hL;
end


end