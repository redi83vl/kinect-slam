%MLESAC Implements modified RANSAC algorithm
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%    [vMask, Model] = MLESAC( mData, ModelFunc, nSampLen, ResidFunc, nIter, dSigma, nEMIter )
%    ----------------------------------------------------------------------------------------------------
%    Arguments:
%           mData - matrix of data, where each column-vector is point
%           ModelFunc - handle to Model Creating function. It must create a
%                   model from nSampLen column-vectors organized in
%                   matrix
%           nSampLen - number of point for ModelFunc
%           ResidFunc - handle to Residuum calculating function. As
%                   argument this function takes model, calculated by
%                   ModelFunc, and matrix of data (all or maybe part of it)
%           nIter - number of iterations for MLESAC algorithm
%           dSigma - sigma value in normal distribution
%           nEMIter - number of iterations in EM algorithm, for estimstion
%                   of mixing parameter
%    Return:
%           vMask - 1s set for inliers, and 0s for outliers
%           Model - approximate model for this data


function [vMask, Model] = MLESAC( mData, ModelFunc, nSampLen, ResidFunc, nIter, dSigma, nEMIter )

% Cheking arguments
if length(size(mData)) ~=2
    error('Data must be organized in column-vecotors massive');
end

nDataLen = size(mData, 2);

% Initialization
Model = NaN;
vMask = zeros([1 nDataLen]);
GlobMaxResid = sqrt(sum( (max(mData') - min(mData')) .^ 2 )); %largest diagonal in this space

dMinPenalty = Inf;

% Main cycle
for i = 1:nIter
    
    % 1. Sampling
    SampleMask = zeros([1 nDataLen]);
    
    % Takes nSampleLen different points
    while sum( SampleMask ) ~= nSampLen
%       SampleMask(randint(1, nSampLen - sum(SampleMask), [1, nDataLen])) = 1;
        ind = ceil(nDataLen .* rand(1, nSampLen - sum(SampleMask)));  
        SampleMask(ind) = 1;
    end    
    Sample = find( SampleMask );
    
    % 2. Creating model
    ModelSet = feval(ModelFunc, mData(:, Sample));
    
    for iModel = 1:size(ModelSet, 3)
        
        CurModel = ModelSet(:, :, iModel);
        
        % 3. Model estimation
        CurResid = abs(feval(ResidFunc, CurModel, mData));
    
        % find mixing parameter, using EM algorithm
        dMix = 0.5; % initialisation
    
        for j = 1:nEMIter
            dResidInlierProb = dMix * exp( -CurResid .^2 / (2 * dSigma ^ 2) ) / (dSigma * sqrt(2 * pi));
            dResidOutlierProb = (1 - dMix) / GlobMaxResid;
            dInlierProb = dResidInlierProb ./ ( dResidInlierProb + dResidOutlierProb );
    
            dMix = mean(dInlierProb);
        end
    
        % find loglkehood of the model
        dResidInlierProb = dMix * exp( -CurResid .^2 / (2 * dSigma ^ 2) ) / (dSigma * sqrt(2 * pi));
        dResidOutlierProb = (1 - dMix) / GlobMaxResid;
        dCurPenalty =  - sum( log( dResidInlierProb + dResidOutlierProb ) );
    
        % 4. The best is selected
        if dMinPenalty > dCurPenalty

            % Save some parameters
            dMinPenalty = dCurPenalty;
            %vMask = dResidInlierProb > dResidOutlierProb;
            vMask = CurResid < 2 * dSigma;
            Model = CurModel;
        end
    end
end

return; 
%END of MLESAC