
% Hoyer's sparsity
% Article: http://www.jmlr.org/papers/v5/hoyer04a.html

% This function 
%   sparseness(A) = 
%               √n (Σ|a_i|/√Σa^2_i) 
%             −−−−−−−−−−−−−−−−−−−−−−−
%                      √n -1
%

function sp = sparsity(A)
%% Returns sparsities of a matrix.
N = size(A,2); % gets number of column
sp = (sqrt(N)-(sum(abs(A'))./sqrt(sum(A'.^2))))/(sqrt(N)-1);

end