function T = Temperature(Q)

global rhoS rhoL CS CL Tc hL;

if (Q < 0)
    % ice 
    T = (Q + rhoS*CS*Tc)/(rhoS*CS);
elseif (Q > rhoL*hL)
    % water
    T = (Q - rhoL*hL+rhoL*CL*Tc)/(rhoL*CL);
else
    T = Tc;
end

end