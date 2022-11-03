function[opt_B,opt_W,  opt_nIter,opt_objhistory, opt_alpha,o,opt_ls]= PenaltyC(trainY, trainX, k, A,B, W, alphamin,alphamax)

if ~exist('B','var')
    B = [];
    W = [];
end
% if isempty(testY) || isempty(testX)
%     testY=[];
%     testX=[];
% end
opt_alpha=alphamin;
options = [];
options.maxIter = 100;
options.alpha=opt_alpha;
options.optimizeB=1;
[opt_B,opt_W, opt_nIter, opt_objhistory,opt_ls]= CVGLGM(trainY, trainX, k, A, options, B, W,10);
% [opt_all,opt_ls]=CalculateObj(testY, opt_B0, testX,opt_B,opt_W,opt_alpha*A,1);
o=opt_ls;
dif=(alphamax-alphamin)/10;
for i=alphamin+0.1:dif:alphamax
    options.alpha=i;
    
%     for r=1:10
%     N0=randi([1 N],N/2)
%     for i=1:N/2
%     n0=randi([length(Y{N0(i)})/2 (length(Y{N0(i)})-1)],1)
%     in=randperm(length(Y{N0(i)}));
%     t=sort(in(1:n0));
%     ct=setdiff([1:length(Y{N0(i)})],t);
%     trainY{N0(i)}=Y{N0(i)}(t);
%     trainX{N0(i)}=X{N0(i)}(:,t);
%     testY{i}=Y{i}(ct);
%     testX{i}=X{i}(:,ct);  
%     end
    [B_final,W_final, nIter_final_1, objhistory_final_1,obj_ls]=CVGLGM(trainY, trainX, k, A, options, B, W,3);
%     [obj_all,obj_ls]=CalculateObj(testY,B0_final, testX,B_final,W_final,i*A,1);
    o=[o,obj_ls];
    if(obj_ls<opt_ls)
        opt_W=W_final;
        opt_B=B_final;
        opt_objhistory=objhistory_final_1;
        opt_nIter=nIter_final_1;
        opt_alpha=i;
        opt_ls=obj_ls;
    end
end
end