function sp = sparsity(A)
%% Returns sparsities of a matrix.
    N = size(A,2); % gets number of column
    sp = (sqrt(N)-(sum(abs(A'))./sqrt(sum(A'.^2))))/(sqrt(N)-1);
end