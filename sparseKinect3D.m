% Input:
% 	im_pts: 2xnpts matrix containing points in the image
% 	W: MxNx3 depth matrix from Kinect device. The generic point (i,j) in
% 	the image has the following coordinates in 3D X=W(i,j,1) Y=W(i,j,2)
% 	and Z=W(i,j,3)
% Output:
% 	w_pts: 3xnpts matrix containing the 3D coordinates of the points from
% 	im_pts. im_pts(k,:) is located in the world in w_pts(k,:)

function [x,y,z] = sparseKinect3D(im_pts,W)
	x=[];
	y=[];
	z=[];
	for k=1:size(im_pts,2)
		
		% 2D
		r=im_pts(1,k);
		c=im_pts(2,k);
		
		% 3D
		x(k)=W(r,c,1);
		y(k)=W(r,c,2);
		z(k)=W(r,c,3);
		
	end
	
	if nargout==1
		x=[x ; y ; z ];
	end
end