function [A,V0]=StochasticBlockModel(NN,P,Diag,Sym)
% function [A,V0] = StochasticBlockModel(NN,P,Diag,Sym)
% Generation of a Stochastic Block Model
%
% Creates a Stochastic Block Model, returns its adjacency matrix.
% This is the classical Stochastic Block Model, with unequal-sized 
% partitions. 
%
% INPUT:
% NN:    vector of  community boundaries
% P:    Matrix of edge probability
% Diag:  if Diag=1, use self-loops; if Diag=0, don't use self-loops
% Sym: if Sym=1, A is sym, if Sym=0 A is not sym
%
% OUTPUT:
% A  adajcency matrix (N-by-N)
% V0 classification vector (N-by-1)
% W  permutation matrix (N-by-N) to put the nodes in order
%
%EXAMPLE
% [A,V0]=StochasticBlockModel([5 10 20 30],[.1 .2 0 .1; ],0,1);
%
K=length(NN);
N=NN(K);NN=[0 NN];

A0=eye(N);

for k1=1:K
    for k2=1:K
    N1=NN(k1)+1;
    N2=NN(k1+1);
    M1=NN(k2)+1;
    M2=NN(k2+1);
    A0(N1:N2,M1:M2)=P(k1,k2);
    end
end

A=eye(N);
if Sym==1
for n1=1:N
    for n2=n1+1:N
        if rand(1)<A0(n1,n2); A(n1,n2)=1; A(n2,n1)=1; end
        
    end
end
else
 for n1=1:N
    for n2=1:N
        if rand(1)<A0(n1,n2); A(n1,n2)=1; end
        
    end
end   
end
%A=A-eye(N);



for k=1:K
	V0(NN(k)+1:NN(k+1),1)=k;
	
end

if Diag~=0
	for n=1:N; A(n,n)=1; end
else	
	for n=1:N; A(n,n)=0; end
end