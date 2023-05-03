function Q = Energy(T)

global Tc hs rhoS rhoL CS CL;

if (T < Tc)
    Q = rhoS*CS*(T-Tc);
else
    Q = rhoL*CL*(T-Tc) + rhoL*hs;
end


end