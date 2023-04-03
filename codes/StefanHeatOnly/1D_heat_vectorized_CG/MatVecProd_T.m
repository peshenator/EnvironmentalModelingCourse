function AT = MatVecProd_T(T,Km,Kp,dQ12)

AT = SpaceDiff_T(T,Km,Kp) + dQ12.*T;

end