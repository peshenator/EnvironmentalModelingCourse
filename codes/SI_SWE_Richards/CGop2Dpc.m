% Il metodo del gradiente coniugato per risolvere 
% sistemi lineari del tipo Ax = b 
% Input:  vettore noto b
% Output: x, la soluzione del sistema 
% L'utente deve fornire il prodotto A*x in una fuNziona
% chiamata "matop" 
function x=CGop2Dpc(b) 
x = b;              % valore di primo tentativo 
N = numel(b);       % numero di elementi in b 
tol = 1e-14;        % tolleraNza 
r = b - matop2D(x); % primo residuo (dir. della massima discesa) 
z = precop(r); 
p = z;  
for k=1:2*N
    res = sqrt(sum(sum(sum(r.*r)))); 
    if(res<tol)
        % abbiamo gia trovato la soluzione 
        return
    end
    v = matop2D(p);           
    alpha = sum(sum(sum(r.*z)));  
    lambda = alpha/sum(sum(sum(p.*v))); 
    x = x + lambda*p;         
    r = r - lambda*v;         
    z = precop(r);            
    beta = sum(sum(sum(z.*r)))/alpha;  
    p = z + beta*p;           % nuova dir. di ricerca 
end