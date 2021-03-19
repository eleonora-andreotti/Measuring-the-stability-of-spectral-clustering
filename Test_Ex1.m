% Example 1: Graph with 17 vertices from the paper. We compute the spectral 
% gaps and the delta_k for k=lb,...,ub (below) and compare them.

lb = 3;
ub = 8;

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

g=specGap(L); % Compute vector with spectral gaps
delta = compute_delta_k(W,lb,ub); % Compute vector with the delta_k

plot(lb:ub,g(lb:ub),'x');
hold on
plot(lb:ub,delta,'o');
legend('Spectral gaps','delta k');
title('Comparison of spectral gaps with the structured measure');
xlabel('k');