function [g] = specGap(L)
% Computes a vector g in R^(n-1) containing the spectral gaps of the matrix
% L in R^n², which is supposed to be the Laplacian of a graph.
s = eig(L);
g = diff(s)/sqrt(2);
end