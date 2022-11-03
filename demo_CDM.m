clc
clear all
addpath('/yourpath/GNMF')
load('/yourpath/simulateddata_1.mat');
%import simulated data using model 1.
%B,W are the coefficients in model 1
[p,N]=size(W);
B0=B*W;
%%plot the dense Y
figure;
for i=1:N
plot(trainY{i}','Color',[0.5,0.5,0.5],'LineWidth',1);
hold on;
end
xlabel('Time','LineWidth',2,'FontSize',14);
ylabel('Response','LineWidth',2,'FontSize',14);

    

%%%-----------------------Experement on Dense Sample------------------%%%%
%%------------Modeling learning-----------------%%
%%IGM

clear BIGM dYIGM
for i=1:N
[BIGM(:,i),bint{i},dYIGM{i}]=regress(trainY{i}',[ones(1,size(trainX{i},2));trainX{i}]');
end
% for i=1:N
% rmse_train(i,1)=sqrt(dYIGM{i}'*dYIGM{i}/length(dYIGM{i}));
% end
%bias
dBIGM=B0-BIGM;
sdB(1)=sum(sum(dBIGM))/N;
vdB(1)=sqrt(sum(diag(dBIGM'*dBIGM))/(N*11));

%%CDM
%%define the parameters
options=[];
options.maxIter = 1000;
options.Converge=0;
options.optimizeB=1;
alphmin=0;
alphmax=1;

%initial value
B_i=[mean(BIGM(:,class==1)')',mean(BIGM(:,class==2)')',mean(BIGM(:,class==3)')'];
w_i=zeros(K,N);
w_i(1,class==1)=1;
w_i(2,class==2)=1;
w_i(3,class==3)=1;

%find the optimal number of clusters
[AIC, K_final, B0_CDM, B_CDM, W_CDM, nIter_CDM, objhistory_CDM] = selectK( trainY', trainX',options, 1, 10);

[B_CDM, W_CDM, nIter_CDM, objhistory_CDM] = CDM (trainY', trainX', K_final, [],options,[],[]);

%goodness of fit
BCDM=B_CDM*W_CDM;
dBCDM=B_1real*W_1real-BCDM;
sdB(2)=sum(sum(dBCDM))/N;
vdB(2)=sqrt(sum(diag(dBCDM'*dBCDM))/(N*11));


%%MEM
Am=ones(N);
%chose the optimal penalty
[B_MEM, W_MEM, nIter_MEM,objhistory_MEM, lambda_MEM,o_MEM]= PenaltyC(trainY', trainX',K_final, Am,[], [], 0.01,10);
options.alpha=lambda_MEM;
[B_MEM, W_MEM, nIter_MEM,objhistory_MEM]= CDM(trainY', trainX',K_final, Am,options,[],[]);
BMEM=B_MEM*W_MEM;
dBMEM=B0-BMEM;
sdB(3)=sum(sum(dBMEM))/N;
vdB(3)=sqrt(sum(diag(dBMEM'*dBMEM))/(N*11));

%%NCDM
Z=W';
%chose the optimal parameters
[B_NCDM, W_NCDM,nIter_NCDM, objhistory_NCDM, objNCDM, lambda, optm, pathobj]= chooseA(trainY,trainX,W',K_final,mean(o_MEM),1,1,20,1);
[B_NCDM2, W_NCDM2,nIter_NCDM2, objhistory_NCDM2, objNCDM2, lambda2, optm2, pathobj2]= chooseA(trainY,trainX,W',K_final,mean(pathobj),1,1,20,3,optm);
[B_NCDM1, W_NCDM1,nIter_NCDM1, objhistory_NCDM1, objNCDM1, lambda1, optm1, pathobj1]= chooseA(trainY,trainX,W',K_final,mean(pathobj),1,1,20,2);
options.alpha=lambda2;
[B0_NCDM, B_NCDM, W_NCDM, nIter_NCDM,objhistory_NCDM, lambda_NCDM,o_NCDM]= PenaltyC(trainY', trainX',K_final, Ak,[], [], 0.01,10);
optm=10;
%%estimate the network
Ak=zeros(N);
%%heat map kernel
id=knnsearch(Z,Z,'k',10);
for i=1:N
    for j=1:optm2
        Ak(i,id(i,j))=exp(-(Z(i,:)-Z(id(i,j),:))*(Z(i,:)-Z(id(i,j),:))'/optm);
    end
end
Ak=(Ak+Ak')/2;
%%0-1
for i=1:N
    Ak(i,id(i,:))=1;
end
Ak=Ak*Ak';
Ak=Ak~=0;

%fit the model
[B_NCDM, W_NCDM, nIter_NCDM,objhistory_NCDM]= CDM(trainY', trainX',K_final, Ak,options,[],[]);
BNCDM=B_NCDM*W_NCDM;
dBNCDM=B0-BNCDM;
sdB(4)=sum(sum(dBNCDM))/N;
vdB(4)=sqrt(sum(diag(dBNCDM'*dBNCDM))/(N*11));

for i=1:N
for j=i:N
Ak(i,j)=exp(-(Z(i,:)-Z(j,:))*(Z(i,:)-Z(j,:))'/optm);
Ak(j,i)=Ak(i,j);
end
end
[objNCDM,lseNCDM,dYNCDM]=CalculateObj(trainY', B0_NCDMk, trainX', B_NCDM, W_NCDM,lambda*Ak);%goodness of fit
rmse_train(:,4)=sqrt(diag(dYNCDM*dYNCDM')/n);
dBNCDM=abs(B0-BNCDM);
sdB(4)=sum(sum(dBNCDM))/nnz(B0);
vdB(4)=mean(var(dBNCDM'));
BNCDM=[B0_NCDM;B_NCDMk*W_NCDM];
dBNCDM=abs(B0-BNCDM);
sdB(4)=sum(sum(dBNCDM))/nnz(B0);
vdB(4)=mean(var(dBNCDM'));
RMSE_train=mean(rmse_train);
figure;
boxplot(rmse_train,'labels',{'IGM','CDM','MEM','NCDM'});
for i=1:N
setot(i,1)=(trainY{i}-mean(trainY{i}))*(trainY{i}-mean(trainY{i}))';
R2(i,1)=1-dYIGM(i,:)*dYIGM(i,:)'/setot(i,1);
R2(i,2)=1-dYCDM(i,:)*dYCDM(i,:)'/setot(i,1);
R2(i,3)=1-dYMEM(i,:)*dYMEM(i,:)'/setot(i,1);
R2(i,4)=1-dYNCDM(i,:)*dYNCDM(i,:)'/setot(i,1);
end
R=mean(R2);

%%------------------------Prediciton------------------%%
nt=5;
for i=1:N
yIGM(i,:)=BIGM(:,i)'*testX{i};
yCDM(i,:)=BCDM(:,i)'*testX{i};
yMEM(i,:)=BMEM(:,i)'*testX{i};
yNCDM(i,:)=BNCDM(:,i)'*testX{i};
end
for i=1:N
rIGM(i,:)=BIGM(:,i)'*testX{i}-testY{i};
rCDM(i,:)=BCDM(:,i)'*testX{i}-testY{i};
rMEM(i,:)=BMEM(:,i)'*testX{i}-testY{i};
rNCDM(i,:)=BNCDM(:,i)'*testX{i}-testY{i};
end
clear rmse_test
rmse_test(:,1)=sqrt(diag(rIGM*rIGM')/nt);
rmse_test(:,2)=sqrt(diag(rCDM*rCDM')/nt);
rmse_test(:,3)=sqrt(diag(rMEM*rMEM')/nt);
rmse_test(:,4)=sqrt(diag(rNCDM*rNCDM')/nt);
figure;
boxplot(rmse_test,'labels',{'IGM','CDM','MEM','NCDM'});

clear RMSE_test
for i=1:5
RMSE_test(i,1)=sqrt(mean(rIGM(:,i).^2))';
RMSE_test(i,2)=sqrt(mean(rCDM(:,i).^2))';
RMSE_test(i,3)=sqrt(mean(rMEM(:,i).^2))';
RMSE_test(i,4)=sqrt(mean(rNCDM(:,i).^2))';
end
RMSE_test=RMSE_test/N;
for i=1:5
    nMSE_test(i,1)=sum(rIGM(:,i).^2)/var(Y(:,20+i));
    nMSE_test(i,2)=sum(rCDM(:,i).^2)/var(Y(:,20+i));
    nMSE_test(i,3)=sum(rMEM(:,i).^2)/var(Y(:,20+i));
    nMSE_test(i,4)=sum(rNCDM(:,i).^2)/var(Y(:,20+i));
end
nMSE_test=sum(nMSE_test)/(5*N);
wR=zeros(1,4);
for i=1:5
wR(1)=corr(Y(:,20+i),yIGM(:,i))*N+wR(1);
wR(2)=corr(Y(:,20+i),yCDM(:,i))*N+wR(2);
wR(3)=corr(Y(:,20+i),yMEM(:,i))*N+wR(3);
wR(4)=corr(Y(:,20+i),yNCDM(:,i))*N+wR(4);
end
wR=wR./(5*N);


%%%%---------------Experiment on Sparse Sample-----------------%%
%%-------------------sparse sampling----------------%%
for i=1:N
    trainXs{i}=trainX{i}(:,1:3:size(trainX{i},2));
    trainYs{i}=trainY{i}(:,1:3:size(trainY{i},2));
end
%%--------------------model learning----------------%%
%%IGM
for i=1:N
[BIGMs(:,i),bints{i},dYIGMs{i}]=regress(trainYs{i}',[ones(1,size(trainXs{i},2));trainXs{i}]');
end
%bias
dBIGMs=B0-BIGMs;
sdBs(1)=sum(sum(dBIGMs))/N;
vdBs(1)=sqrt(sum(diag(dBIGMs*dBIGMs'))/(N*11));
%CDM
options.alpha = 0;
[AICs, K_finals, B0_CDMs, B_CDMs, W_CDMs, nIter_CDMs, objhistory_CDMs] = selectK( trainYs', trainXs',options, 1, 10);
[B_CDMs, W_CDMs, nIter_CDMs, objhistory_CDMs] = CDM (trainYs', trainXs', K, [],options,[],[]);
%goodness of fit
BCDMs=B_CDMs*W_CDMs;
dBCDMs=B0-BCDMs;
sdBs(2)=sum(sum(dBCDMs))/N;
vdBs(2)=sqrt(sum(diag(dBCDMs*dBCDMs'))/(N*11));
%%MEM
[B_MEMs, W_MEMs, nIter_MEMs,objhistory_MEMs, lambda_MEMs,o_MEMs]= PenaltyC(trainYs', trainXs',K_finals, Am,[], [], 0.01,10);
options.alpha=lambda_MEMs;
options.maxIter = 1000;
[B_MEMs, W_MEMs, nIter_MEMs,objhistory_MEMs]= CDM(trainYs', trainXs',K_finals, Am,options,[],[]);
BMEMs=B_MEMs*W_MEMs;
dBMEMs=B0-BMEMs;
sdBs(3)=sum(sum(dBMEMs))/N;
vdBs(3)=sqrt(sum(diag(dBMEMs*dBMEMs'))/(N*11));
%%NCDM
Z=W';
optm=1.1;
for i=1:N
for j=i:N
Ak(i,j)=exp(-(Z(i,:)-Z(j,:))*(Z(i,:)-Z(j,:))'/optm);
Ak(j,i)=Ak(i,j);
end
end
[B_NCDMs, W_NCDMs,nIter_NCDMs, objhistory_NCDMs, objNCDMs, lambdas, optms, pathobjs]= chooseA(trainYs,trainXs,W',K_finals,min(o_MEMs),1,1,20,1);
[B_NCDMs2, W_NCDMs2,nIter_NCDMs2, objhistory_NCDMs2, objNCDMs2, lambdas2, optms2, pathobjs2]= chooseA(trainYs,trainXs,W',K_finals,min(o_MEMs),1,1,20,3,optms);

[B_NCDMs, W_NCDMs, nIter_NCDMs,objhistory_NCDMs, lambda_NCDMs,o_NCDMs]= PenaltyC(trainYs', trainXs',K, Ak,[], [], 0.01,10);
id=knnsearch(Z,Z,'k',optms2);
Aks=zeros(N);
for i=1:N
    for j=1:optms2
        Aks(i,id(i,j))=exp(-(Z(i,:)-Z(id(i,j),:))*(Z(i,:)-Z(id(i,j),:))'/optms);
    end
end
Aks=(Aks+Aks')/2;

    Aks=zeros(N);
    for i=1:N
    Aks(i,id(i,:))=1;
    end
    
options.alpha=lambdas2;
[B_NCDMs, W_NCDMs, nIter_NCDMs,objhistory_NCDMs]= CDM(trainYs', trainXs',K_finals, Aks,options,[],[]);
BNCDMs=B_NCDMs*W_NCDMs;
dBNCDMs=B0-BNCDMs;
sdBs(4)=sum(sum(dBNCDMs))/N;
vdBs(4)=sqrt(sum(diag(dBNCDMs*dBNCDMs'))/(N*11));
nt=5;
for i=1:N
yIGMs(i,:)=BIGMs(:,i)'*testX{i};
yCDMs(i,:)=BCDMs(:,i)'*testX{i};
yMEMs(i,:)=BMEMs(:,i)'*testX{i};
yNCDMs(i,:)=BNCDMs(:,i)'*testX{i};
end
for i=1:N
rIGMs(i,:)=BIGMs(:,i)'*testX{i}-testY{i};
rCDMs(i,:)=BCDMs(:,i)'*testX{i}-testY{i};
rMEMs(i,:)=BMEMs(:,i)'*testX{i}-testY{i};
rNCDMs(i,:)=BNCDMs(:,i)'*testX{i}-testY{i};
end
rmse_tests(:,1)=sqrt(diag(rIGMs*rIGMs')/nt);
rmse_tests(:,2)=sqrt(diag(rCDMs*rCDMs')/nt);
rmse_tests(:,3)=sqrt(diag(rMEMs*rMEMs')/nt);
rmse_tests(:,4)=sqrt(diag(rNCDMs*rNCDMs')/nt);
figure;
boxplot(rmse_tests,'labels',{'IGM','CDM','MEM','NCDM'});


%%%----------------------Prediction-------------------%%%
clear RMSE_tests
for i=1:5
RMSE_tests(i,1)=sqrt(mean(rIGMs(:,i).^2))';
RMSE_tests(i,2)=sqrt(mean(rCDMs(:,i).^2))';
RMSE_tests(i,3)=sqrt(mean(rMEMs(:,i).^2))';
RMSE_tests(i,4)=sqrt(mean(rNCDMs(:,i).^2))';
end
RMSE_test=RMSE_test/N;
for i=1:5
    nMSE_tests(i,1)=sum(rIGMs(:,i).^2)/var(Y(:,20+i));
    nMSE_tests(i,2)=sum(rCDMs(:,i).^2)/var(Y(:,20+i));
    nMSE_tests(i,3)=sum(rMEMs(:,i).^2)/var(Y(:,20+i));
    nMSE_tests(i,4)=sum(rNCDMs(:,i).^2)/var(Y(:,20+i));
end
nMSE_tests=sum(nMSE_tests)/(5*N);
wRs=zeros(1,4);
for i=1:5
wRs(1)=corr(Y(:,20+i),yIGMs(:,i))*N+wRs(1);
wRs(2)=corr(Y(:,20+i),yCDMs(:,i))*N+wRs(2);
wRs(3)=corr(Y(:,20+i),yMEMs(:,i))*N+wRs(3);
wRs(4)=corr(Y(:,20+i),yNCDMs(:,i))*N+wRs(4);
end
wRs=wRs./(5*N);