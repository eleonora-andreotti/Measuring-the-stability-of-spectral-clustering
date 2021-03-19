function [eps_opt,eps_lb,eps_ub,E_opt] = newtonBisection(W,k,eps_lb,eps_ub,eps0,m_max,tol)
% Input parameters:
% W :       Weight matrix of the initial graph
% k :       Index of delta_k to be computed
% eps0 :    Starting value of epsilon: An educated guess would be the k-th
%           spectral gap since delta_k is larger
% eps_lb:   Initial lower bound of eps_opt
% eps_ub:   Initial upper bound of eps_opt
% tol:      Tolerance value
% m_max:    Maximum number of Newton-bisection iterations

gamma = 0.8; % damping parameter, possibly as input argument!

%% Check input parameters
if nargin<7 % tolerance not specified
  tol = 1e-7;
  if nargin<6 % maximum iteration not specified
    m_max = 1000;
    if nargin < 5 % starting value of epsilon not specified, use eps_lb
      if nargin>3
        eps0 = eps_lb;
      end
      if nargin<3
        eps0 = 0;
        eps_lb = 0;
        eps_ub = norm(Lap(W),'fro');
      if nargin<2
        error('Not enough input arguments');
      end
      end
    end
  end
end

%% Initializations
[V,~] = eig(Lap(W));    % Compute eigenvectors and eigenvalues of Lap(W)
%V(:,1) = ones(size(V,1),1);
%V = orth(V(:,1:k+1));
x = V(:,k+1);           % ... to compute G_eps(0)
y = V(:,k);
G_eps = LapStar(x*x'-y*y',W>0);
E0 = -G_eps/norm(Lap(G_eps),'fro');

m = 0; % Initialize iteration counter
[E,f,f_der] = innerIteration(W,k,eps0,E0); % Perform first iteraion

%% Start method
eps_new = eps0;
while m<=m_max
    eps_old = eps_new;
    if abs(f)<tol % eigenvalues coalesce => eps>=eps*, perform bisection step
        eps_ub = min(eps_ub,eps_old);
        eps_new = (eps_lb+eps_ub)/2;
    else % eigenvalues don't coalesce <=> f_eps>0 => perform newton step
        eps_lb = max(eps_lb,eps_old);
        eps_new = eps_old-gamma*f/f_der;
        %E0 = eps_old/eps_new*E;
        E0 = E;
    end
    
    % check if Newton step yields a new epsilon which excesses the bounds
    if eps_new<eps_lb || eps_new>eps_ub
        eps_new = (eps_lb+eps_ub)/2;
    end
    % Stop if maximum iteration number is reached or the remaining
    % intervall is smaller then the tolerance
    if m==m_max || eps_ub-eps_lb < tol
        if m==m_max
            fprintf('Maximum number of iterations reached.\n');
        end
        eps_opt = eps_new;
        m = m_max+1;
    else
        m = m+1;
    end
    [E,f,f_der] = innerIteration(W,k,eps_new,E0);
end
E_opt = E;
end