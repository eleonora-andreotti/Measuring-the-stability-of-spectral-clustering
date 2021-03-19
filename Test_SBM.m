clear; close; clc
% Example of a SBM with 4 blocks of equal size. A block of size n is
% represented by a 2-vertex graph whose nodes are connected by weights
% connected in the vector alpha.
% The random perturbation between the first two blocks is represented by
% small parameters mu.
mu = 35:0.5:55;
% For small mu, we expect a clustering into 4 clusters to be more robust
% than in 3 clusters. We compute the thresholds mu*, where the spectral gap
% and the new measure choose a 3-clustering to be more robust.
N = 30;
alpha = 100*ones(1,N);
W = createBlockGraph(alpha);
% W(3:4,1:2) = mu*eye(2);
% W(1:2,3:4) = mu*eye(2); % represents random perturbation
k_opt_g = zeros(length(mu),1);
k_opt_d = k_opt_g;

for k=1:length(mu)
    W(3:4,1:2) = mu(k)*eye(2);
    W(1:2,3:4) = mu(k)*eye(2); % represents random perturbation
    L = Lap(W);
    g = specGap(L);
    g=g(1:2*N-2);
    [~,k_opt_g(k)] = max(g);
    delta = zeros(N,1);
    delta(N-1:N) = compute_delta_k(W,N-1,N);
    [~,k_opt_d(k)] = max(delta);
end
plot(mu,k_opt_g,'x');
hold on;
plot(mu,k_opt_d,'o');
title('Optimal number of clusters for different perturbation parameters');
legend('k_{opt} for spectral gaps','k_{opt} for \delta_k');
xlabel('\mu');
axis([mu(1),mu(end),N-2,N+1])