%MSAC Implements modified RANSAC algorithm
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%    [vMask, Model] = MSAC( mData, ModelFunc, nSampLen, ResidFunc, nIter, dThreshold, iAdaptive )
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
%           nIter - number of iterations for MSAC algorithm
%           dThreshold - threshold for residuum
%           iAdaptive  - 0 if not adaptive, >=1 if adaptive, iAdaptive *
%           number of iterations for current 
%    Return:
%           vMask - 1s set for inliers, and 0s for outliers
%           Model - approximate model for this data


function [vMask, Model, nResIter] = MSAC( mData, ModelFunc, nSampLen, ResidFunc, nIter, dThreshold, iAdaptive )

if nargin < 7 iAdaptive = 0; end;
if nargin < 8 FitParam = []; end;

nResIter = nIter;
    
% Cheking arguments
if length(size(mData)) ~=2
    error('Data must be organized in column-vecotors massive');
end

nDataLen = size(mData, 2);

if( nDataLen < nSampLen )
    error('Not enough data to compute model function');
end

% Initialization
Model = NaN;
vMask = zeros([1 nDataLen]);

dMinPenalty = Inf;


% Main cycle
for i = 1:nIter
    
    % 1. Sampling
    SampleMask = zeros([1 nDataLen]);
    
    % Takes nSampleLen different points
    while sum( SampleMask ) ~= nSampLen
%        SampleMask(randint(1, nSampLen - sum(SampleMask), [1, nDataLen])) = 1;
        ind = ceil(nDataLen .* rand(1, nSampLen - sum(SampleMask)));  
        SampleMask(ind) = 1;
    end    
    Sample = find( SampleMask );
    
    % 2. Creating model
    ModelSet = feval(ModelFunc, mData(:, Sample));

    for iModel = 1:size(ModelSet, 3)
        
        CurModel = ModelSet(:, :, iModel);
    
        % 3. Model estimation    
        CurResid    = abs(feval(ResidFunc, CurModel, mData));
        dCurPenalty = sum(min(CurResid, dThreshold));
    
        % 4. The best is selected
        if dMinPenalty > dCurPenalty

            % Save some parameters
            dMinPenalty = dCurPenalty;
            vMask = (CurResid < dThreshold);
            Model = CurModel;
        end
    end
    
    %adaptive finish. Calculate probability of success based on current
    %iteration and outlier level
    if iAdaptive ~= 0 
        dOutliers = (size( vMask,2 ) - sum( vMask )) / size( vMask,2 );
        p = sac_success_prob( nSampLen, dOutliers, round(i / iAdaptive) );
        if p > 0.95
            nResIter = i;
            break;
        end;
    end;
    
end

return; 
%END of MSAC