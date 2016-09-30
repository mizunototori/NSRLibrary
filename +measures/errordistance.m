% Frobenius norm as error distance between A and B

function error = errordistance(A,B)
error = norm((A-B),'fro');
end
