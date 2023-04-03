function res = ResidualOuter(T,Km,Kp,rhs)

 res = Q(T) + SpaceDiff_T(T,Km,Kp) - rhs;
%  res = 0*Q(T) + SpaceDiff_T(T,Km,Kp) - 0*rhs;
 
end