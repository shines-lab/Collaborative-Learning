function [B0_final, B_final,W_final, nIter_final, objhistory_final,o,obj_ls0]= CVGLGM(Y, X, k, A, options, B, W,n)
N=length(Y);
obj_ls0=10^6;
o=[];
for r=1:n
    clear trainY trainX testY testX
%     N0=randperm(N*10);
%     N0=sort(N0(1:N));
%     index=zeros(N,10);
%     for i=1:length(N0)
%         n1=mod(N0(i),10);
%         n2=ceil(N0(i)/10);
%         if(n1==0)
%             n1=10;
%         end
%         index(n2,n1)=1;
%     end
% %     n0=randi([length(Y{N0(i)})/2 (length(Y{N0(i)})-1)],1)
% %     in=randperm(length(Y{N0(i)}));
% %     t=sort(in(1:n0));
% %     ct=setdiff([1:length(Y{n2})],n1);
%     for i=1:N
%     trainY{i}=Y{i}(index(i,:)==0);
%     trainX{i}=X{i}(:,index(i,:)==0);
%     testY{i}=Y{i}(index(i,:)==1);
%     testX{i}=X{i}(:,index(i,:)==1);  
for i=1:N
    in=randperm(length(Y{i})-1)+1;
    trainY{i}=Y{i}([1,sort(in(1:end-1))]);
    trainX{i}=X{i}(:,[1,sort(in(1:end-1))]);
    testY{i}=Y{i}(in(end));
    testX{i}=X{i}(:,in(end));    
end
    [B,W, nIter, objhistory]= CDM(trainY', trainX', k, A, options, B, W);
    [obj_all,obj_ls]=CalculateObj(testY',[], testX',B,W,options.alpha*A,1);
    obj_ls=obj_ls/N;
    o=[o,obj_ls];
    if(obj_ls<obj_ls0)
        B_final=B;
        W_final=W;
        nIter_final=nIter;
        objhistory_final=objhistory;
        obj_ls0=obj_ls;
    end
end
o=mean(o);