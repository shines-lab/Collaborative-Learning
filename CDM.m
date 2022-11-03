function [B_final, W_final, nIter_final, objhistory_final] = CDM (Y, T, k, A, options, B, W)
% Collaborative Degradation Model (NCDM)
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

if ~isfield(options,'error')
    options.error = 1e-5;
end
if ~isfield(options, 'maxIter')
    options.maxIter = [];
end

if ~isfield(options,'nRepeat')
    options.nRepeat = 10;
end

if ~isfield(options,'minIter')
    options.minIter = 30;
end

if ~isfield(options,'meanFitRatio')
    options.meanFitRatio = 0.1;
end

if ~isfield(options,'alpha')
    if(isempty(A))
        options.alpha = 0;
    else
        options.alpha=100;
    end
end

if ~isfield(options,'optimizeB')
    options.optimizeB = 0;
end

nSmp = size(Y,1);

if isfield(options,'alpha_nSmp') && options.alpha_nSmp
    options.alpha = options.alpha*nSmp;    
end


if ~isfield(options,'Optimization')
    options.Optimization = 'Multiplicative';
end

if ~exist('B','var')
    B = [];
    W = [];
end

switch lower(options.Optimization)
    case {lower('Multiplicative')} 
         [B_final, W_final, nIter_final, objhistory_final] = CDM_Multi(Y, T, k, A, options, B, W);
    otherwise
        error('optimization method does not exist!');
end


    
        