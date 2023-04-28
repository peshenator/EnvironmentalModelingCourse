% Thomas algorithm for the solution of tridiagonal systems 
function x=Thomas(a,b,c,d) 
N=numel(d);     % compute the number of elements in the system 
c(1)=c(1)/b(1); 
d(1)=d(1)/b(1); 
% Part I: forward elimination 
for i=2:N
    gamma=1/(b(i)-c(i-1)*a(i)); 
    c(i)=c(i)*gamma;
    d(i)=(d(i)-a(i)*d(i-1))*gamma; 
end
% Part II: back substitution 
x = zeros(1,N); % allocate row vector 
x(N) = d(N);    
for i=N-1:-1:1 % reverse loop 
  x(i) = d(i)-c(i)*x(i+1);   
end
