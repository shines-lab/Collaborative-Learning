function [B0_NCDM, B_NCDM, W_NCDM,nIter_NCDM, objhistory_NCDM, objall_w, lambda, optm, pathobj]= chooseA(trainY,trainX,Z,k,objall_w,mmin,mdif,mmax,scheme,h)
options=[];
options.maxIter = 100;
% Bi=B_CDM;
% Wi=W_CDM;
alphamin=0;
alphamax=1;
N=length(trainY);
pathobj=objall_w;
 if(scheme==1)
for m=mmin:mdif:mmax
for i=1:N
    for j=i:N
        A(i,j)=exp(-(Z(i,:)-Z(j,:))*(Z(i,:)-Z(j,:))'/m);
        A(j,i)=A(i,j);
    end
end
       
% for i=1:N
%     in=randperm(length(Y{i})-1)+1;
%     trainY{i}=Y{i}([1,sort(in(1:end-1))]);
%     trainX{i}=X{i}(:,[1,sort(in(1:end-1))]);
%     testY{i}=Y{i}(in(end));
%     testX{i}=X{i}(:,in(end));    
% end
[opt_B0, opt_B, opt_W, opt_nIter,opt_objhistory, opt_alpha,o, objall_test]= PenaltyC(trainY', trainX',k, A,[], [], alphamin,alphamax);
% [objall_test,objCDM_test]=CalculateObj(testY',opt_B0, testX',opt_B,opt_W,opt_alpha*A,1);
% objCDM_test=sqrt(objCDM_test./(N*p));
if(objall_test<objall_w)
    B0_NCDM=opt_B0;
    B_NCDM=opt_B;
    W_NCDM=opt_W;
%     Inter_NCDM=Inter_test;
    objhistory_NCDM=opt_objhistory;
    nIter_NCDM=opt_nIter;
    objall_w=objall_test
%     objCDM_w=objCDM_test
    lambda=opt_alpha
    optm=m
%     Bi=B_NCDM;
%     Wi=W_NCDM;
end
pathobj=[pathobj, objall_test];
end
 elseif(scheme==2 | scheme==3)
 for m=mmin:mdif:mmax
id=knnsearch(Z,Z,'k',m);
A=zeros(N);
if(scheme==2)
for i=1:N
    A(i,id(i,:))=1;
end
A=A+A';
A=A~=0;
else
for i=1:N
    for j=1:m
        A(i,id(i,j))=exp(-(Z(i,:)-Z(id(i,j),:))*(Z(i,:)-Z(id(i,j),:))'/h);
    end
end
A=(A+A')/2;
end
       
% for i=1:N
%     in=randperm(length(Y{i})-1)+1;
%     trainY{i}=Y{i}([1,sort(in(1:end-1))]);
%     trainX{i}=X{i}(:,[1,sort(in(1:end-1))]);
%     testY{i}=Y{i}(in(end));
%     testX{i}=X{i}(:,in(end));    
% end
[opt_B0, opt_B, opt_W, opt_nIter,opt_objhistory, opt_alpha,o, objall_test]= PenaltyC(trainY', trainX',k, A,[], [], alphamin,alphamax);
% [objall_test,objCDM_test]=CalculateObj(testY',opt_B0, testX',opt_B,opt_W,opt_alpha*A,1);
% objCDM_test=sqrt(objCDM_test./(N*p));
if(objall_test<objall_w)
    B0_NCDM=opt_B0;
    B_NCDM=opt_B;
    W_NCDM=opt_W;
%     Inter_NCDM=Inter_test;
    objhistory_NCDM=opt_objhistory;
    nIter_NCDM=opt_nIter;
    objall_w=objall_test
%     objCDM_w=objCDM_test
    lambda=opt_alpha
    optm=m
%     Bi=B_NCDM;
%     Wi=W_NCDM;
end
pathobj=[pathobj, objall_test];
 end   
 end
end