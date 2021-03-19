function [E,F_eps,F_eps_der] = innerIteration(W,k,epsilon,E0,tol)
% Computes the optimal perturbation matrix E such that
% F_eps(E) = lambda(k+1)(Lap(W+epsilonE))-lambda(k)(Lap(W+epsilonE))
% is minimized under the constraint that Lap(E) has unit Frobenius norm.
% 
% returns also the derivative of F_eps(E*) with respect to epsilon, i. e.
% F_eps_der = d/d eps [ F(eps,E*(eps)) ]
if nargin<5
    tol = 1e-9;
    if nargin<4 % Initial value not specified
        [V,~] = eig(Lap(W));    % Compute eigenvectors and eigenvalues of Lap(W)
        V(:,1) = ones(size(V,1),1);
        V = orth(V(:,1:k+1));
        x = V(:,k+1);           % ... to compute G_eps(0)
        y = V(:,k);
        G_eps = LapStar(x*x'-y*y',W>0);
        E0 = -G_eps/norm(Lap(G_eps),'fro');
    end
end

P_W = W>0;              % matrix containing ones where W has edges
h = 0.1;                % initial step size
h_min = tol/100;
h_max = h;
E = E0;

%% perform first iteration
L = Lap(W+epsilon*E0);
[V,lambda] = eig(L);
lambda = diag(lambda);
x = V(:,k+1);
y = V(:,k);
F_eps = lambda(k+1)-lambda(k);
F_eps_old = F_eps+1;

ITER = 1; % iteration counter, can be deleted

while abs(F_eps-F_eps_old)>max(abs(tol*h*F_eps),tol/1000) && abs(F_eps)>tol
    % stopping criteria could be optimized, very restrictive!
    
    % Compute functional G_eps and E_dot
    G_eps = LapStar(x*x'-y*y',P_W);
    LLE = LapStar(Lap(E),P_W);
    kappa = sp_frob(G_eps,LLE)/sp_frob(LLE,LLE);
    E_dot = -G_eps+kappa*LLE;
  %  E_dot = E_dot/norm(E_dot,'fro'); % normalize E_dot
    E_old = E;
    E = E+h*E_dot;
    
    % Normalize E1 such that Lap(E1) has norm 1
    E = E/norm(Lap(E),'fro');
    
    % Compute new eigenvalues
    L = Lap(W+epsilon*E);
    [V,lambda] = eig(L);
    lambda = diag(lambda);
    x = V(:,k+1);
    y = V(:,k);
    F_eps_old = F_eps;
    F_eps = lambda(k+1)-lambda(k);
    h_inc = 1;
    ITER = ITER+1;
    while F_eps > F_eps_old && h > h_min
        h_inc = 0;
        h = h/2;
        E = E_old+h*E_dot;
        E = E/norm(Lap(E),'fro');
        L = Lap(W+epsilon*E);
        [V,lambda] = eig(L);
        lambda = diag(lambda);
        x = V(:,k+1);
        y = V(:,k);
        F_eps = lambda(k+1)-lambda(k);
        ITER = ITER+1;
    end
    if h_inc == 1 && 1.2*h<h_max
        h = 1.2*h;
    end
end
G_eps = P_W.*Sym(LapStar(x*x'-y*y',P_W));
LLE = LapStar(Lap(E),P_W);
F_eps_der = -norm(G_eps,'fro')/norm(LLE,'fro');
end