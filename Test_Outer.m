% Test outer iteration

% Example 1: Graph with 17 vertices from the paper. We compute the spectral 
% gaps and the delta_k for k=lb,...,ub (below) and compare them.

lb = 2;
ub = 5;

W = [0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
     0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 1;
     0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
     0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 1;
     0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
     0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 1;
     0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
     0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 1;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
     0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
W = W+W';
L = Lap(W);

for k=lb:ub
    [eps_opt,eps_lb,eps_ub] = newtonBisection(W,k)
end

eps_opts = compute_delta_k(W)