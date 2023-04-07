function res = ResidualOuter(T,Kmx,Kpx,Kmy,Kpy,rhs)

res = Q(T) + SpaceDiff_T(T,Kmx,Kpx,Kmy,Kpy) - rhs;
 
end