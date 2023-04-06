% Jacobian dQ2/dT of Q2(T)
% because dQ = dQ1 - dQ2 --> dQ2 = dQ1 - dQ
function jac = dQ2(T)

jac = dQ1(T) - dQ(T);

end
