%LOMSAC Implements MSAC, modified with LO procedure
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%    function [vMask, Model] = LOMSAC( mData, ModelFunc, nSampLen, ResidFunc, nIter, dThreshold, LOFunc, nLOSampLen )
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
%           LOFunc - function of local optimization of model 
%           nLOSampLen - number of samples to supply to LOFunc
%    Return:
%           vMask - 1s set for inliers, and 0s for outliers
%           Model - approximate model for this data


function [vMask, Model, nResIter] = LOMSAC( mData, ModelFunc, nSampLen, ResidFunc, nIter, dThreshold, LOModelFunc, nLOSampLen, iAdaptive )

if nargin < 9
    iAdaptive = 0;
end;
nResIter = nIter;

% Cheking arguments
if length(size(mData)) ~=2
    error('Data must be organized in column-vecotors massive');
end

nDataLen = size(mData, 2);

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
        SampleMask(randint(1, nSampLen - sum(SampleMask), [1, nDataLen])) = 1;
    end    
    Sample = find( SampleMask );
    
    % 2. Creating model
    ModelSet = feval(ModelFunc, mData(:, Sample));

    for iModel = 1:size(ModelSet, 3)
        
        CurModel = ModelSet(:, :, iModel);

    
        % 3. Model estimation    
        CurResid = abs(feval(ResidFunc, CurModel, mData));
        dCurPenalty = sum(min(CurResid, dThreshold));
        InlierMask = CurResid < dThreshold;
        
        
    
        % 4. The best is selected
        if dMinPenalty > dCurPenalty
            
            if sum( InlierMask ) <= (nLOSampLen * 4) 
              %if too few of inliers, don't to LO step  
              dMinPenalty = dCurPenalty;
              vMask       = InlierMask;
              Model       = CurModel;
            else
              % LO step when number of inliers is big enough  
              while sum( SampleMask & InlierMask ) ~= nLOSampLen
                   SampleMask( randint( 1, nLOSampLen - sum( SampleMask & InlierMask ), [1, nDataLen])) = 1;
              end
            
              Sample     = find( SampleMask & InlierMask );
              LOModel    = feval( LOModelFunc, mData( :, Sample ) );
              LOResid    = abs( feval( ResidFunc, LOModel, mData ) );
              dLOPanalty = sum( min( LOResid, dThreshold ) );

              % Save some parameters
              dMinPenalty = dLOPanalty;
              vMask       = (LOResid < dThreshold);
              Model       = LOModel;
            end
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
%END of LOMSAC