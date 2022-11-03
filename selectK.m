function [AIC, K_final, B_final, W_final, nIter_final, objhistory_final] = selectK( Y, X,options,Kmin, Kmax)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
AIC=[];
oldAIC=10^10;
for k=Kmin:1:Kmax
    [B_CDM, W_CDM, nIter_CDM, objhistory_CDM] = CDM (Y, X, k, [],options,[],[]);
    [objall,objCDM,dYCDM]=CalculateObj(Y,[],X,B_CDM,W_CDM,[],1);
    newAIC=objall+2*2*k;
    if(newAIC<oldAIC)
        B_final=B_CDM;
        W_final=W_CDM;
        nIter_final=nIter_CDM;
        objhistory_final=objhistory_CDM;
        K_final=k;
        oldAIC=newAIC;
    end
    AIC=[AIC,newAIC];
    
end
end

