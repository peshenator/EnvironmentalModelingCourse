function AT = MatVecProd_T(T,Kmx,Kpx,Kmy,Kpy,dQ12)

AT = SpaceDiff_T(T,Kmx,Kpx,Kmy,Kpy) + dQ12.*T;

end