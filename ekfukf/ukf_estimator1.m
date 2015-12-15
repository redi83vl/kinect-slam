function [x,P]=ukf_estimator1(x,P,f,nx,y,h,ny)
	iter=2;
	for k=1:iter
% 		prediction
		[x,P] = UKF_PREDICT1(x,P,f,nx);
%		correction
		[x,P] = UKF_UPDATE1(x,P,y,h,ny);
	end
		
end