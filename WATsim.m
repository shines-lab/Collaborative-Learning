function [W,A,class]=WAsim(N,K,C,sigmaK,p)
%1.Generate K pieces of Gaussian Distribution
%2.Generate W from the mixture of the K pieces
%3.Generate A from the correlation matrix of W
%maxT: the maximul time point
%N: namber of samples
%n: Nx1 vector contains the number of time points of each sample
%K: number of pieces
%C: number of classes
%sigmaK: Kx1 cells, each cell contains the covariance matrix of each
%distribution piece
%p: the probability of mixture pieces
if(size(sigmaK,1)~=K || length(p)~=K)
    error('check the dimension of sigmaK');
end
for i=1:N
%     s=randperm(maxT);
%     t=sort(s(1:n(i)));
%     T{i}(1,:)=ones(1,n(i));
%     T{i}(2,:)=t;
%     T{i}(3,:)=t.^2;
    class(i)=discreteinvrnd(p,1,1);
    if(size(sigmaK{class(i)},1)~=C)
        error('check the dimension of each sigma')
    else
        W(:,i)=abs(mvnrnd(zeros(C,1),sigmaK{class(i)}));
        W(:,i)=W(:,i)./sum(W(:,i));
    end
end
A=normc(W)'*normc(W);


     
      
        
        
        
        
        
        
        
        
function X = discreteinvrnd(p,m,n)

X = zeros(m,n); % Preallocate memory
for i = 1:m*n
    u = rand;
    I = find(u < cumsum(p));
    X(i) = min(I);
end

    
    