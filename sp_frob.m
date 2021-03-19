function [y] = sp_frob(A,B)
% Computes the inner product of matrices A and B corresponding to the
% frobenius norm. We assume that A and B are square matrices of the same
% size.
y = sum(sum(A.*B));
end

