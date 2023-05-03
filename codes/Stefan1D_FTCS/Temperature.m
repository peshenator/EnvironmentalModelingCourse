function T = Temperature(Q)

global rhoS rhoL CS CL Tc hs;

if (Q < 0)
    % ice 
    T = (Q + rhoS*CS*Tc)/(rhoS*CS);
elseif (Q > rhoL*hs)
    % water
    T = (Q - rhoL*hs+rhoL*CL*Tc)/(rhoL*CL);
else
    T = Tc;
end

end