function [T,inl]=rigidTransformation(X,Y,method)
	M=[X;Y];
	
	if nargin==2
		T=simpleRigidTransformation(M);
		inl = ones(1,size(M,2));
	elseif strcmp(method,'ransac')
		[T,inl]=ransacRigidTransformation(M);
		T=simpleRigidTransformation(M(:,inl));
	elseif strcmp(method,'rransac')
		[T,inl]=rransacRigidTransformation(M);
		T=simpleRigidTransformation(M(:,inl));
	elseif strcmp(method,'lmeds')
		T=lmedsRigidTransformation(M);
		inl = ones(1,size(M,2));
	else
		error('??? Unkmown method %s',method);
	end
	
end

% Input:
%     X,Y: two set of points represented as 3-by-K matrices. It is assumed
%     that the point X(:,k) in the original frame correspondes to the the
%     point Y(:,k) in the destination frame.
% Output:
%     T: the rigid transfromation that connects the two sets, i.e.:
%     X(:,j)=T*Y(:,j).

function T=simpleRigidTransformation(M)
	X=M(1:3,:);
	Y=M(4:6,:);
	
	T=absor(X,Y);
	T=inv(T.M);
end

function [T,inl]=ransacRigidTransformation(M)
	mData = M;
	ModelFunc = @simpleRigidTransformation;
	nSampLen = 3;
	ResidFunc = @errorRes;
	nIter = 10;
	dThreshold = .1;
	[vMask, Model] = RANSAC( mData, ModelFunc, nSampLen, ResidFunc, nIter, dThreshold );% ransac from GML library
	
	T = Model;
	inl = vMask;
end

function [T,inl]=rransacRigidTransformation(M)
	mData = M;
	ModelFunc = @simpleRigidTransformation;
	nSampLen = 3;
	ResidFunc = @errorRes;
	nIter = 10;
	dThreshold = .1;
	nTestLen = size(M,2)/2;
	[vMask, Model] = RRANSAC( mData, ModelFunc, nSampLen, ResidFunc, nIter, dThreshold, nTestLen );% rransac from GML library
	T = Model;
	inl = vMask;
end

function T=lmedsRigidTransformation(M)
	mData = M;
	ModelFunc = @simpleRigidTransformation;
	nSampLen = 3;
	ResidFunc = @errorRes;
	nIter = 10;
	Model = LMEDS( mData, ModelFunc, nSampLen, ResidFunc, nIter );
	T = Model;
end


function err = errorRes(T,M)
	X = [ M(1:3,:) ; ones(1,size(M,2)) ];
	Y = [ M(4:6,:) ; ones(1,size(M,2)) ];
% 	residual Errccor
	err=sum( T*Y - X , 1 ).^2;
end













