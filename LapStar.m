function [L] = LapStar(V,P_E)
% Computes the adjoint Laplacian of a matrix V, where P_E is the sparsity
% pattern (an indicator matrix with 0 and 1).
n = size(V,1);
X = diag(V)*ones(1,n)-V;
L = P_E.*Sym(X);
end

