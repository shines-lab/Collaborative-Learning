clc
clear all
%%%-----------------model 1----------------------%%%
K=3;
sigK=1000;
N=100;
sigY=5;
p=[0.4,0.3,0.3];
% B=[400,0,-0.5;400,-40,1;400,0,-1;400,0,0;400,-20,0]';
%%Simulate the s
Bs=[10,-0.01;5,-0.1;2,-0.5]';
for i=1:K
s=1*ones(K,1);
s(i)=s(i)+sigK;
sigmaK{i}=spdiags(s,0,K,K);
end
[W,A,class]=WATsim(N,K,K,sigmaK',p);
clear Tall
for i=1:N
Tall{i}=[ones(1,25);1:25];
end
S=MCsim(W,Tall,Bs,0.1);
figure;
plot(S','Color',[0.5,0.5,0.5]);
hold on;
h1=plot([1:25],mean(S(class==1,:)),'r',[1:25],mean(S(class==2,:)),'b',[1:25],mean(S(class==3,:)),'g','LineWidth',2);
xlabel('Time','LineWidth',2);
ylabel('Response','LineWidth',2);
legend(h1,'Class 1','Class 2','Class 3');
%%simulate cov(X)
sigmaX=spdiags(10*ones(4,1),0,4,4);
BX(:,2)=normrnd(0,50,10,1);
BX(:,3)=normrnd(0.5,0.01,10,1);
BX(:,1)=5*binornd(2,0.5,10,1)+2;
BX(:,4)=5*binornd(2,0.5,10,1)+2;
clear fun
f=(0:0.01:1);
    for m=1:10
        for t=1:length(f)
            fun(m,t)=BX(m,1)*(1+exp(-BX(m,2)*(f(t)-BX(m,3))))^-1+BX(m,4);
        end
    end
    figure;
plot(fun','LineWidth',2);
xlim=([0 100]);
set(gca,'XTickLabel',{'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'})
legend('Marker1','Marker2','Marker3','Marker4','Marker5','Marker6','Marker7','Marker8','Marker9','Marker10');
%%simulate X
nS=normr(S);
clear X
for i=1:N
    for m=1:10
        for t=1:25
            e=normrnd(0,2,1);
            X{i}(m,t)=BX(m,1)*(1+exp(-BX(m,2)*(nS(i,t)-BX(m,3))))^-1+BX(m,4)+e;
        end
    end
end
%%simulate Y
clear B
 sigmaB1=spdiags([100,100,100,100,100,0.1,0.1,0.1,0.1,0.1]',0,10,10);
 sigmaB2=spdiags(ones(1,10)',0,10,10);
 sigmaB3=spdiags([0.1,0.1,0.1,0.1,0.1,100,100,100,100,100]',0,10,10);
B(:,1)=mvnrnd(zeros(10,1)',sigmaB1,1)';
B(:,2)=mvnrnd(-5*ones(10,1)',sigmaB2,1)';
B(:,3)=mvnrnd(zeros(10,1)',sigmaB3,1)';
% %%way 2
% for i=1:10
%     sigmaB=spdiags(randi(5,1)*ones(3,1),0,3,3);
%     B(i,:)=mvnrnd(zeros(3,1)',sigmaB,1)';
% end

Y=MCsim(W,X,B,sigY);
endT=20;
clear B0;
for i=1:N
    B0(1,i)=-min(Y(i,:));    
    Y(i,:)=Y(i,:)-min(Y(i,:));
end
B0=[B0;B*W];
clear trainX trainY testX testY validX validY
for i=1:N
trainX{i}=X{i}(:,1:endT);
% validX{i}=X{i}(2:3,ct(i,:));
trainY{i}=Y(i,1:endT);
% validY{i}=Y(i,ct(i,:));
testX{i}=[ones(1,5);X{i}(:,21:25)];
testY{i}=Y(i,21:25);
end
N0=randperm(N);
N0=sort(N0(1:N/2));
for i=1:length(N0)
    n=randi([5 10],1);
    s=randperm(endT);
t=sort(s(1:n));
ct=setdiff([1:20],t);
trainX{N0(i)}=X{N0(i)}(:,t);
validX{N0(i)}=X{N0(i)}(:,ct);
trainY{N0(i)}=Y(N0(i),t);
validY{N0(i)}=Y(N0(i),ct);  
end
%%plot the dense Y
figure;
for i=1:N
plot(trainY{i}','Color',[0.5,0.5,0.5],'LineWidth',1);
hold on;
end
% h1=plot([1:20],mean(Y(class==1,1:20)),'r',[1:20],mean(Y(class==2,:)),'b',[1:25],mean(Y(class==3,:)),'g','LineWidth',2);
xlabel('Time','LineWidth',2,'FontSize',14);
ylabel('Response','LineWidth',2,'FontSize',14);
% legend(h1,'Class 1','Class 2','Class 3');

    
%%%-------------------model 2--------------------%%%
K=5;
sigK=1000;
N=100;
sigY=5;
p=[0.2,0.2,0.2,0.2,0.2];
% B=[400,0,-0.5;400,-40,1;400,0,-1;400,0,0;400,-20,0]';
%%Simulate the s
Bs=[10,-0.01;10,-0.2;7,-0.5;2,-0.5;2,-0.01]';
for i=1:K
s=1*ones(K,1);
s(i)=s(i)+sigK;
sigmaK{i}=spdiags(s,0,K,K);
end
[W,A,class]=WATsim(N,K,K,sigmaK',p);
clear Tall
for i=1:N
Tall{i}=[ones(1,25);1:25];
end
S=MCsim(W,Tall,Bs,0.1);
figure;
plot(S','Color',[0.5,0.5,0.5]);
hold on;
h1=plot([1:25],mean(S(class==1,:)),'r',[1:25],mean(S(class==2,:)),'b',[1:25],mean(S(class==3,:)),'g','LineWidth',2);
xlabel('Time','LineWidth',2);
ylabel('Response','LineWidth',2);
legend(h1,'Class 1','Class 2','Class 3');
%%simulate cov(X)
sigmaX=spdiags(10*ones(4,1),0,4,4);
BX(:,2)=normrnd(0,50,10,1);
BX(:,3)=normrnd(0.5,0.01,10,1);
BX(:,1)=5*binornd(2,0.5,10,1)+2;
BX(:,4)=5*binornd(2,0.5,10,1)+2;
clear fun
f=(0:0.01:1);
    for m=1:10
        for t=1:length(f)
            fun(m,t)=BX(m,1)*(1+exp(-BX(m,2)*(f(t)-BX(m,3))))^-1+BX(m,4);
        end
    end
    figure;
plot(fun','LineWidth',2);
xlim=([0 100]);
set(gca,'XTickLabel',{'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'})
legend('Marker1','Marker2','Marker3','Marker4','Marker5','Marker6','Marker7','Marker8','Marker9','Marker10');
%%simulate X
nS=normr(S);
clear X
for i=1:N
    for m=1:10
        for t=1:25
            e=normrnd(0,2,1);
            X{i}(m,t)=BX(m,1)*(1+exp(-BX(m,2)*(nS(i,t)-BX(m,3))))^-1+BX(m,4)+e;
        end
    end
end
%%simulate Y
clear B
 sigmaB1=spdiags([100,100,100,100,100,0.01,0.01,0.01,0.01,0.01]',0,10,10);
 sigmaB2=spdiags(ones(1,10)',0,10,10);
 sigmaB3=spdiags([0.01,0.01,0.01,0.01,0.01,100,100,100,100,100]',0,10,10);
 sigmaB4=spdiags([0.1,100,0.1,100,0.1,100,0.1,100,0.1,100]',0,10,10);
 sigmaB5=spdiags(0.1*ones(1,10)',0,10,10);
B(:,1)=mvnrnd(zeros(10,1)',sigmaB1,1)';
B(:,2)=mvnrnd(-5*ones(10,1)',sigmaB2,1)';
B(:,3)=mvnrnd(zeros(10,1)',sigmaB3,1)';
B(:,4)=mvnrnd(5*ones(10,1)',sigmaB4,1)';
B(:,5)=mvnrnd(-1*ones(10,1)',sigmaB5,1)';
B(1:5,5)=0;
B(6:10,4)=0;
% %%way 2
% for i=1:10
%     sigmaB=spdiags(randi(5,1)*ones(3,1),0,3,3);
%     B(i,:)=mvnrnd(zeros(3,1)',sigmaB,1)';
% end

Y2=MCsim(W,X,B,sigY);
endT=20;
clear B0;
for i=1:N
   % B0(1,i)=max(0,-min(Y(i,:)));    
    Y2(i,:)=Y2(i,:)-min(Y2(i,:));
end
B0=[B0;B*W];
clear trainX trainY testX testY validX validY
for i=1:N
trainX{i}=X{i}(:,1:endT);
% validX{i}=X{i}(2:3,ct(i,:));
trainY{i}=Y(i,1:endT);
% validY{i}=Y(i,ct(i,:));
testX{i}=[ones(1,5);X{i}(:,21:25)];
testY{i}=Y2(i,21:25);
end
N0=randperm(N);
N0=sort(N0(1:N/2));
for i=1:length(N0)
    n=randi([5 10],1);
    s=randperm(endT);
t=sort(s(1:n));
ct=setdiff([1:20],t);
trainX{N0(i)}=X{N0(i)}(:,t);
validX{N0(i)}=X{N0(i)}(:,ct);
trainY{N0(i)}=Y(N0(i),t);
validY{N0(i)}=Y(N0(i),ct);  
end
%%plot the dense Y
figure;
for i=1:N
plot(trainY{i}','Color',[0.5,0.5,0.5],'LineWidth',1);
hold on;
end
% h1=plot([1:20],mean(Y(class==1,1:20)),'r',[1:20],mean(Y(class==2,:)),'b',[1:25],mean(Y(class==3,:)),'g','LineWidth',2);
xlabel('Time','LineWidth',2,'FontSize',14);
ylabel('Response','LineWidth',2,'FontSize',14);
% legend(h1,'Class 1','Class 2','Class 3');
