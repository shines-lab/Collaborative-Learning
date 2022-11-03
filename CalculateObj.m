function [obj, obj_NMF,dY] = CalculateObj(Y, B0,T, B, W, A,index)
if ~exist('index')
    index=0;
end
    MAXARRAY = 500*1024*1024/8; % 500M. You can modify this number based on your machine's computational power.
if ~exist('B0')
     B0 = zeros(1,size(Y,1));
end
%     if ~exist('dVordU','var')
%         dVordU = 1;
%     end
    dV = [];
    mFea = size(Y,1);
    mn = numel(A);
    nBlock = ceil(mn/MAXARRAY);

    if mn < MAXARRAY
        obj_NMF=0;
        for i=1:mFea
            if(isempty(Y{i})==0)
                if(index==1)
                dY{i} = (Y{i}-B0(i))'-T{i}'*B*W(:,i);
                obj_NMF = obj_NMF+dY{i}'*dY{i};
                else
                dY(i,:) = Y{i}-B0(i)-(T{i}'*B*W(:,i))';
                obj_NMF = obj_NMF+dY(i,:)*dY(i,:)';
                end
            end
        end
%         if deltaVU
%             if dVordU
%                 dV = dX'*U + L*W;
%             else
%                 dV = dX*V;
%             end
%         end
    else
        obj_NMF = 0;
%         if deltaVU
%             if dVordU
%                 dB = zeros(size(B));
%             else
%                 dW = zeros(size(W));
%             end
%         end
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
%             if deltaVU
%                 if dVordU
%                     dV(smpIdx,:) = dX'*U;
%                 else
%                     dV = dU+dX*V(smpIdx,:);
%                 end
%             end
        end
%         if deltaVU
%             if dVordU
%                 dV = dV + L*V;
%             end
%         end
    end
    if isempty(A)
        obj_Lap = 0;
    else
%         obj_Lap =  sum(sum((W*L).*W));
obj_Lap=0;
for i=1:mFea
    for j=1:mFea
        if(j~=i)
            obj_Lap=obj_Lap+(W(:,i)-W(:,j))'*(W(:,i)-W(:,j))*A(i,j);
        end
    end
end

    end
    obj = obj_NMF+obj_Lap;
    

