function [A]=simA (K,N,sigma)

x1 = 0.3+0.4*rand;
x2n0 = rand(1,K-1);
sumx2n0 = sum(x2n0);
index = N*[ x1, (1-x1)*x2n0/sumx2n0];
index=round(index);
index(K)=N-sum(index(1:K-1));
A=[];
for i=1:K
    cov{i}=normrnd(0,sigma,index(i));
    cov{i}=abs(corr(cov{i}));
    A=blkdiag(A,cov{i});
end