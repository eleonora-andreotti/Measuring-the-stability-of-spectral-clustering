function [delta] = compute_delta_k(W,k_min,k_max)
% Computes
%           delta_k = min||L(W)-L(W*)||
% under the condition that the k-th and (k+1)-th eigenvalue of L(W*)
% coalesce and under the constraints that W* is a weight matrix with the
% same sparsity pattern as W

% delta = compute_delta_k(W) computes is equivalent to
% compute_delta_k(W,1,size(W)-1)
% delta = compute_delta_k(W,k) computes delta_k
% delta = compute_delta_k(W,k1,k2) computes a vector delta containing
% delta_k for k between k_min and k_max

if nargin < 3
    if nargin<2
        disp('test');
        k_min=1;
        k_max = size(W,1)-1;
    else
        k_max = k_min;
    end
end

delta = zeros(k_max-k_min+1,1);
gaps = specGap(Lap(W));     % Compute spectral gaps...
for j=1:length(delta)
    fprintf('%d out of %d deltas computed.\n',j,length(delta));
    k = k_min+j-1;
    eps0 = gaps(k);         % ... to use k-th spectral gap as an initial
                            % value for epsilon.
    [eps_opt,~,~,~] = newtonBisection(W,k_min+j-1,eps0,norm(Lap(W),'fro'),eps0);
    delta(j) = eps_opt;
end
end