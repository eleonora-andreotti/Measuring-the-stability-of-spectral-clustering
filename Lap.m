function [L] = Lap(W)
% Computes the Laplacian matrix associated with the weight matrix W.
n = size(W,1);
assert(n==size(W,2),'Input matrix not symmetric.');
L = diag(W*ones(n,1))-W;
end