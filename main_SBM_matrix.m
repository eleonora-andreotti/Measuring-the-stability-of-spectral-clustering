clc
close all
clear


%NN=[30,60,70,120];
NN=[25,50,75,100];
Diag=0;
Sym=1;
pdiag=.9;
pout=.02;
pint=3*.02;

P=[pdiag,pint,pout,pout;
    pint,pdiag,pout,pout;
    pout,pout,pdiag,pint;
    pout,pout,pint,pdiag]
    
[A,V0]=StochasticBlockModel(NN,P,Diag,Sym);

% A=A.*randi(1,dim);A=A+A';
% A=A/norm(A,'fro');