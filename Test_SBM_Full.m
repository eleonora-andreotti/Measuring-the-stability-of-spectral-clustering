clc
close all
clear

N = 6;
f = @(MU0,k)(-MU0/(N-1)*k+MU0+MU0/(N-1));
mu0 = 0.05:0.05:1;
%mu0 = 1
%NN=[30,60,70,120];
%NN=[25,50,75,100];
NN = zeros(1,N);
for k=1:N
    NN(1,k) = 25*k;
end


%% Start routine
pin = ones(N,1);
K = length(mu0); % Number of complete computations to be done
pout = zeros(N-1,1);
k_opt_g = zeros(K,1);
k_opt_d = zeros(K,1);
delta = zeros(N,K);


for k=1:K
    fprintf('Computation at %d percent.\n',100*k/K);
    % Determine parameters for block perturbations
    for l=1:N-1
        pout(l,1) = f(mu0(k),l);
    end
    P = diag(pin)+diag(pout,1)+diag(pout,-1);
    [W,~]=StochasticBlockModel(NN,P,0,1);
    
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

