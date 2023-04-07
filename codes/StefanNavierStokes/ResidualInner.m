function [res,dQ12] = ResidualInner(Talpha1,Talpha,Kmx,Kpx,Kmy,Kpy,rhs)

res = Q1(Talpha1) - Q2(Talpha) - dQ2(Talpha).*(Talpha1 - Talpha) - rhs;
% add M*T
res = res + SpaceDiff_T(Talpha1,Kmx,Kpx,Kmy,Kpy);

dQ12 = dQ1(Talpha1) - dQ2(Talpha);

end