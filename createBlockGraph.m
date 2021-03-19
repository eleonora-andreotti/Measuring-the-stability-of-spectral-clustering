function [W] = createBlockGraph(alpha)
% Creates weight matrix of a block graph used to represent stochastic block
% models.

% For alpha = [a1,...aN], a graph with 2N vertices is created such that it
% has N clusters of 2 vertices each, which are connected with weight aj,
% j=1,.,N
% Example: For alpha = [1 2], a graph with weight matrix
% [0 1 0 0
%  1 0 0 0
%  0 0 0 2
%  0 0 2 0]
% is created.

n = length(alpha);
W = zeros(2*n,2*n);
for k=1:n
    W(2*k-1,2*k) = alpha(k);
    W(2*k,2*k-1) = alpha(k);
end
end