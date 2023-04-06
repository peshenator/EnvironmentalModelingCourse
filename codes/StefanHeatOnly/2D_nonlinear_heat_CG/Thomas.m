% Thomas algorithm for the solution of 
% linear tridiagonal systems 
function x=Thomas(a,b,c,d) 
c(1) = c(1)/b(1); 
d(1) = d(1)/b(1); 
N = numel(d); % get the number of elements in d 
% Part I: forward elimination 
for i=2:N
    g = 1/( b(i)-c(i-1)*a(i) );
    c(i) = c(i)*g; 
    d(i) = (d(i)-a(i)*d(i-1))*g; 
end
% Part II: back substitution 
x = zeros(1,N); % allocate x as a column vector 
x(N) = d(N); 
for i=N-1:-1:1
    x(i) = d(i) - c(i)*x(i+1); 
end




