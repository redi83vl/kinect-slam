function [X,Y,Z,R,G,B]=mxNiPointCloud(KinectHandles,sz)
	
	if nargin==1
		sz=1;
	end
	
	W = mxNiDepthRealWorld(KinectHandles); W = permute(W,[2 1 3]);
	I = mxNiPhoto(KinectHandles); I = permute(I,[3,2,1]);
	
	if nargin==2
		W=imresize(W,sz);
		I=imresize(I,sz);
	end
	
	if nargout==2
		X=W;
		Y=I;
		return;
	end
	
	%3D points coordinates from depth map
	X=W(:,:,1);X=X(:);
	Y=W(:,:,2);Y=Y(:);
	Z=W(:,:,3);Z=Z(:);
	
	%3D points colors from the image
	R=I(:,:,1);R=R(:);
	G=I(:,:,2);G=G(:);
	B=I(:,:,3);B=B(:);
	
	inl = logical( (Z > 0.3) .* (Z<10000) );% good points
	X=X(inl)/1000;	Y=Y(inl)/1000;	Z=Z(inl)/1000;% 3D points coordinates
	R=double(R(inl));	G=double(G(inl));	B=double(B(inl));% 3D points colors
end