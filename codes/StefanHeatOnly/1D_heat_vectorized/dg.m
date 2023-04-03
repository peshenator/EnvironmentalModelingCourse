% approximate derivative of the function g, 
% using a central finite difference 
function dy=dg(gamma)
epsilon=1e-7; 
dy = ( g(gamma+epsilon)-g(gamma-epsilon) )/(2*epsilon); 

