% Thomas algorithm for the solution of tridiagonal systems 
% Input: 
%   a = vector of lower diagonal elements
%   b = vector of diagonal elements
%   c = vector of upper diagonal elements
%   d = known right hand side vector 
% Output: 
%   x = solution vector 
% Function name: Thomas
% File name :    Thomas.m (compulsory) 
function x=Thomas(a,b,c,d)
N = length(d);  % get the number of elements (problem size) 
% Part I: forward elimination 
c(1)=c(1)/b(1); 
d(1)=d(1)/b(1); 
for i=2:N
    gamma=1/(b(i)-c(i-1)*a(i)); 
    c(i) = c(i)*gamma; 
    d(i) = (d(i)-a(i)*d(i-1))*gamma; 
end
% Part II: back substitution 
x(N)=d(N); 
for i=N-1:-1:1
    x(i) = d(i)-c(i)*x(i+1); 
end
x=x(:);  % force MATLAB to produce a column vector 



