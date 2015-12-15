%ZHANGSAC Implements modified RANSAC algorithm (with MSAC cost function)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%    [vMask, Model] = ZHANGSAC( mData, ModelFunc, nSampLen, ResidFunc, nIter, dThreshold, dRadius )
%    ---------------------------------------------------------------------------------
%    Arguments:
%           mData - matrix of data, where each column-vector is point
%           ModelFunc - handle to Model Creating function. It must create a
%                   model from nSampLen column-vectors organized in
%                   matrix
%           nSampLen - number of point for ModelFunc
%           ResidFunc - handle to Residuum calculating function. As
%                   argument this function takes model, calculated by
%                   ModelFunc, and matrix of data (all or maybe part of it)
%           nIter - number of iterations for ZHANGSAC algorithm
%           dThreshold - threshold for residuum
%           dRadius - radius of hypersphere for subsampling
%    Return:
%           vMask - 1s set for inliers, and 0s for outliers
%           Model - approximate model for this data


function [vMask, Model] = ZHANGSAC( mData, ModelFunc, nSampLen, ResidFunc, nIter, dThreshold, dRadius )

% Cheking arguments
if length(size(mData)) ~=2
    error('Data must be organized in column-vecotors massive');
end

nDataLen = size(mData, 2);
dRadius2 = dRadius ^ 2;

% Initialization
Model = NaN;
vMask = zeros([1 nDataLen]);

dMinPenalty = Inf;

% Main cycle
for i = 1:nIter
    
    % 1. Sampling
    Sample = zeros([1 nSampLen]);
    Sample(1) = ceil(nDataLen .* rand(1)); %randint(1, 1, [1, nDataLen]);
    
    % Mask of all points out sphere
    SampleMask = int8( sum((mData - mData(:, Sample(1)) * ones([1 nDataLen])) .^ 2) > dRadius2 );
    SampleMask(Sample(1)) = 1;
    
    % Sample fails if not enouph points
    nSpherePoints = sum(SampleMask);
    if nSpherePoints < nSampLen
        continue;
    end
    
    % Takes nSampleLen different points
    SphereInd = find(SampleMask);
    SampleMask(Sample(1)) = 2;  % Alredy has been selected    
    
    while sum( SampleMask ) ~= (nSampLen + nSpherePoints)
%        ind = randint(1, nSampLen + nSpherePoints - sum(SampleMask), [1, nSpherePoints]);
        ind = ceil(nSpherePoints .* rand(1, nSampLen + nSpherePoints - sum(SampleMask))); 
        SampleMask(SphereInd(ind)) = 2;
    end
    Sample = find( SampleMask == 2 );
    
    % 2. Creating model
    ModelSet = feval(ModelFunc, mData(:, Sample));
    
    for iModel = 1:size(ModelSet, 3)
        
        CurModel = ModelSet(:, :, iModel);
        
        % 3. Model estimation    
        CurResid = abs(feval(ResidFunc, CurModel, mData));
        dCurPenalty = sum(min(CurResid, dThreshold));
    
        % 4. The best is selected
        if dMinPenalty > dCurPenalty

            % Save some parameters
            dMinPenalty = dCurPenalty;
            vMask = (CurResid < dThreshold);
            Model = CurModel;
        end
    end
end

return; 
%END of ZHANGSAC