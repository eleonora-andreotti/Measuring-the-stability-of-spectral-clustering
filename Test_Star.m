clear; close; clc
% Example of a SBM with N blocks of equal size. A block of size n is
% represented by a 2-vertex graph whose nodes are connected by weights
% connected in the vector alpha.

%% Input here
N = 8;  % Number of blocks
%f = @(MU0,k)(MU0/k); % Function for decrease in perturbations
%f = @(MU0,k)(MU0/(1.2^(k-1)));
%f = @(MU0,k)(-MU0/(N-1)*k+MU0+MU0/(N-1));
%f = @(MU0,k)(MU0-MU0/N^2*k^2);
mu0 = 5:5:200;

%% Start routine
K = length(mu0); % Number of complete computations to be done
mus = zeros(N,1);
k_opt_g = zeros(K,1);
k_opt_d = zeros(K,1);
delta = zeros(N,K);


for k=1:K
    fprintf('Computation at %d percent.\n',100*k/K);
    % Determine parameters for block perturbations
    for l=1:N
        mus(l,1) = f(mu0(k),l);
    end
    alpha = 100*ones(1,N);
    W = createBlockGraph(alpha); % Create Graph
    
    % Add perturbartions to graph
    for l=1:N-1
        W(1:2,2*l+1:2*l+2) = mus(l)*eye(2);
        W(2*l+1:2*l+2,1:2) = mus(l)*eye(2);
    end
    
    % Compute spectral gaps and determine maximum
    L = Lap(W);
    g = specGap(L);
    g = g(1:N);
    [~,k_opt_g(k)] = max(g);
    
    % Compute stuctured robustness measure and determine maximum
%    delta(1:2,k) = [0;0];
    delta(1:N,k) = compute_delta_k(W,1,N);
    [~,k_opt_d(k)] = max(delta(:,k));
end

figure(1)
plot(mu0,k_opt_g,'x');
hold on;
plot(mu0,k_opt_d,'o');
title('Optimal number of clusters for different perturbation parameters');
legend('k_{opt} for spectral gaps','k_{opt} for \delta_k');
xlabel('\mu');
axis([mu0(1),mu0(end),1,N+1]);

figure(2)
plot(mu0,delta(N,:));
hold on;
plot(mu0,delta(N-1,:));
plot(mu0,delta(N-2,:));
plot(mu0,delta(2,:));
legend('\delta_8','\delta_7','\delta_6','\delta_2');

save
