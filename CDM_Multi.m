function [B_final, W_final, nIter_final, objhistory_final] = CDM_Multi(Y, T, k, A, options, B, W)
% Collaborative Degradation Modeling (NCDM) with
%          multiplicative update
% 
%  
% Notation:
% Y ... (mFea x 1) cells
%       mFea ... number of subjects
% Y{i}  (1 x nSmp(i))  outcomes 
%       nSmp(i)  ... number of time points of subject i
% T ... (mFea x 1) cells
% T{i}  (p x nSmp(i)) predictor matrix i.e [1,t,t^2]'
% k ... number of latent group
% A ... weight matrix of the network 
%
% options ... Structure holding all settings
%
% You only need to provide the above four inputs.
% 
% Output:
% B_final ... (p x k) parameter matrix (latent cluster)
% W_final ... (k x mFea) weight matrix (membership vector)
%
% References:
% [1] Lin, Y., Liu, K., Byon, E., Qian, X., and Huang, S., ?Domain-Knowledge Driven Cognitive
% Degradation Modeling for Alzheimer`s Disease?, SDM 2015.
%
%
%

%   version 1.0 --Oct/2015 
%
%   Written by Ying Lin (linyeliana DOT ie AT gmail.com)
%

differror = options.error;
maxIter = options.maxIter;
nRepeat = options.nRepeat;
minIter = options.minIter - 1;
optimizeB=options.optimizeB;
if ~isempty(maxIter) && maxIter < minIter
    minIter = maxIter;
end
meanFitRatio = options.meanFitRatio;

alpha = options.alpha;
if isempty(A) && alpha>0
    error('network matrix missing');
end
Norm = 0;

mFea=size(Y,1);
n=size(T{1},1);
U=[];%initial U
objhistory_final=[];
objhistory=[];
for i=1:mFea
    nSmp(i)=size(Y{i},2);
    if(size(T{i},2)==size(Y{i},2))
    U=[U,T{i}*Y{i}'];%used in update W
    else
        error('size of Y and T do not match');
    end
end

if alpha > 0 
    A = alpha*A;
    DCol = full(sum(A,2));
    D = spdiags(DCol,0,mFea,mFea);
    L = D - A;
    if isfield(options,'NormW') && options.NormA
        D_mhalf = spdiags(DCol.^-.5,0,mFea,mFea) ;
        L = D_mhalf*L*D_mhalf;
    end
else
    L = [];
    A= [];
end

selectInit = 1;
if isempty(B)
    B = rand(n,k);
    W = abs(rand(k,mFea));
    objhistory=[];
else
    nRepeat = 1;
end

W = NormalizeW(W, Norm);
if nRepeat == 1
    selectInit = 0;
    minIter = 0;
    if isempty(maxIter)
        objhistory = CalculateObj(Y, T, B, W, L);
        meanFit = objhistory*10;
    else
        objhistory_final = CalculateObj(Y, T, B, W, L);
        if isfield(options,'Converge') && options.Converge
            objhistory = CalculateObj(Y, T, B, W, L);
        end
    end        
else
    if isfield(options,'Converge') && options.Converge
        error('Not implemented!');
    end
end




tryNo = 0;
nIter = 0;

while tryNo < nRepeat   
    tryNo = tryNo+1;
    maxErr = 1;
    while(maxErr > differror || maxErr<0)
        % ===================== update W ========================
        BV=[];
        BU=[];
        for i=1:mFea
            BU=[BU,B'*T{i}*Y{i}'];
            BV=[BV,B'*T{i}*T{i}'*B*W(:,i)];
        end  
        if alpha > 0
            AW = W*A;
            DW = W*D;
            
            BU = BU + AW;
            BV = BV + DW;
        end
        
        W = W.*(BU./max(BV,1e-10));
%         W = NormalizeW(W, Norm);
        % ===================== update B ========================
        %Update using LSE closed form
        if(optimizeB==0)
            XX=zeros(n*k);
            XY=zeros(n*k,1);
            for i=1:mFea
                for j=1:nSmp(i)
                    r=[];
                    for v=1:n
                        r=[r,T{i}(v,j)*W(:,i)'];
                    end%use LSE to get B
                    X{i}(j,:)=r;
                end
                XX=XX+X{i}'*X{i}+2*eye(n*k);
                XY=XY+X{i}'*Y{i}';
            end
            b=pinv(XX)*XY;%b is a 3k x 1 vector
            b=b';
            B=[];
            for v=1:n
                B=[B;b((v-1)*k+1:v*k)];
            end%transfer vector b to matrix B
        %Update using QP
        %You need to install the cvx or other QP package in you matlab
        else
            cvx_begin
                variable B(n,k);
                expression optB(mFea,1);
            for i=1:mFea
                optB(i,1)=norm(Y{i}'-T{i}'*B*W(:,i),2);
            end
            minimize(sum(optB))
            subject to
                for i=1:mFea
                T{i}'*B>=0;
                end
            cvx_end
        end
        % ===================== update objhistory ========================
        nIter = nIter + 1;
        if nIter > minIter
            if selectInit
                objhistory = CalculateObj(Y, T, B, W, L);
                maxErr = 0;
            else
             %o=[o,CalculateObj(Y, T, B, W, L)];
                if isempty(maxIter)
                    newobj = CalculateObj(Y, T, B, W, L);
                    objhistory = [objhistory newobj]; %#ok<AGROW>
                    meanFit = meanFitRatio*meanFit + (1-meanFitRatio)*newobj;
                    maxErr = (meanFit-newobj)/meanFit;
                else
                    if isfield(options,'Converge')
                        newobj = CalculateObj(Y, T, B, W, L);
                        objhistory = [objhistory newobj]; %#ok<AGROW>
                   end
                    maxErr = 1;
                    if nIter >= maxIter
                        maxErr = 0;
                        if isfield(options,'Converge') && options.Converge
                        else
                            objhistory = [objhistory 0];
                        end
                    end
                end
            end
        end
    end
    
    if tryNo == 1 && nRepeat~=1
        B_final = B;
        W_final = W;
        nIter_final = nIter;
        objhistory_final = objhistory;
    else
       if objhistory(end) < objhistory_final(end)
           B_final = B;
           W_final = W;
           nIter_final = nIter;
           objhistory_final = objhistory;
       end
    end

    if selectInit
        if tryNo < nRepeat
            %re-start
            B = rand(n,k);
            W = abs(rand(k,mFea));
            
            W = NormalizeW(W, Norm);
            nIter = 0;
        else
            tryNo = tryNo - 1;
            nIter = minIter+1;
            selectInit = 0;
            B = B_final;
            W = W_final;
            objhistory = objhistory_final;
            meanFit = objhistory*10;
        end
    end
end



%==========================================================================

function [obj, dV] = CalculateObj(Y, T, B, W, A)
    MAXARRAY = 500*1024*1024/8; % 500M. You can modify this number based on your machine's computational power.
    dV = [];
    mFea = size(Y,1);
    mn = numel(A);
    nBlock = ceil(mn/MAXARRAY);

    if mn < MAXARRAY
        obj_NMF=0;
        for i=1:mFea
            dY{i} = Y{i}'-T{i}'*B*W(:,i);
            obj_NMF = obj_NMF+dY{i}'*dY{i};
        end
    else
        obj_NMF = 0;
        PatchSize = ceil(nSmp/nBlock);
        for i = 1:nBlock
            if i*PatchSize > mFea
                smpIdx = (i-1)*PatchSize+1:mFea;
            else
                smpIdx = (i-1)*PatchSize+1:i*PatchSize;
            end
            for j=1:length(smpIdx)
                dY{j} = Y{smpIdx(j)}'-T{smpIdx(j)}'*B*W(:,smpIdx(j));
                obj_NMF = obj_NMF+dY{j}'*dY{j};
            end
        end
    end
    if isempty(A)
        obj_Lap = 0;
    else
        obj_Lap=0;
        Norm=1;
        nW = NormalizeW(W, Norm);
        obj_Lap=sum(sum((nW*A).*nW));
    end
    obj = obj_NMF+obj_Lap;
    




function [W] = NormalizeW(W, Norm)
    N = size(W,2);
    if Norm == 2
            norms = max(1e-15,sqrt(sum(W.^2,1)))';
            W = W*spdiags(norms.^-1,0,N,N);
    else
            norms = max(1e-15,sum(abs(W),1))';
            W = W*spdiags(norms.^-1,0,N,N);

    end

        