% Planted partitioning model
clear; close; clc;
alpha = [20 20 20];
p = 0.1:0.05:0.4;
n = length(p);
g2 = zeros(n,1);
g3 = zeros(n,1);
d2 = zeros(n,1);
d3 = zeros(n,1);
for k=1:n
    k
    W = createPPMGraphFirstGroup(alpha,p(k));
    g = specGap(Lap(W));
    spy(W)
    drawnow;
    g2(k,1)=g(2);
    g3(k,1)=g(3);
    delta = compute_delta_k(W,2,3);
    d2(k,1)=delta(1);
    d3(k,1)=delta(2);
end

subplot(1,2,1)
plot(p,g2,'Marker','o');
hold on
plot(p,g3,'Marker','o');
title('Spectral gaps for different values of p');
legend('g_2','g_3');
subplot(1,2,2)
plot(p,d2,'Marker','x');
hold on
plot(p,d3,'Marker','o');
title('Structured measure for different values of p');
legend('d_2','d_3')

profile viewer