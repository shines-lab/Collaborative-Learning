function[Y] = MCsim(W,T,B,sigmaE)
%T has N cells; T{i} is 3*ni matrix with 1st row is 1, 2nd row is time,
%3rd row is time^2;
[N,K]=size(W');
if(length(T)~=N || size(B,2)~=K)
    error('check the size of INPUT');
end
Y=[];
for i=1:N
    ni=size(T{i},2);
    e=normrnd(0,sigmaE(i),ni,1);
    Y=[Y;W(:,i)'*B'*T{i}+e'];
end
