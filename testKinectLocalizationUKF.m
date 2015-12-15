%%
clc;
clear all;
close all;

%% Start the Kinect Process
xml_config_file='C:\Users\Redjan\Desktop\workspace\imports\matlab\win_kin_mex\Config\SamplesConfig.xml';
context = mxNiCreateContext(xml_config_file);
mxNiSetProperty(context, 'MinDepthValue', 10, 'MinDepthValue', 5000);
K = [ 526.37013657, 0.00000000, 313.68782938
	0.00000000, 526.37013657, 259.01834898 
	0.00000000, 0.00000000, 1.00000000 ];
%% Start LIBVISIO Matcher
% matching parameters
nms_n                  = 5;   % non-max-suppression: min. distance between maxima (in pixels)
nms_tau                = 50;  % non-max-suppression: interest point peakiness threshold
match_binsize          = 50;  % matching bin width/height (affects efficiency only)
match_radius           = 200; % matching radius (du/dv in pixels)
match_disp_tolerance   = 1;   % du tolerance for stereo matches (in pixels)
outlier_disp_tolerance = 5;   % outlier removal: disparity tolerance (in pixels)
outlier_flow_tolerance = 5;   % outlier removal: flow tolerance (in pixels)
multi_stage            = 1;   % 0=disabled,1=multistage matching (denser and faster)
half_resolution        = 0;   % 0=disabled,1=match at half resolution, refine at full resolution
refinement             = 0;   % refinement (0=none,1=pixel,2=subpixel)

matcherMex('init',nms_n,nms_tau,match_binsize,match_radius,match_disp_tolerance,outlier_disp_tolerance,outlier_flow_tolerance,multi_stage,half_resolution,refinement);

%%

deltaT=.1;

first = 2;
step = 1;
last = 1000;
index = first;

T=eye(4);

tic;
[W0,I0]=mxNiPointCloud(context);
I0 = rgb2gray(I0);
matcherMex('push',I0);


subplot(3,3,[4 5 7 8]), hold on, grid on, axis equal, view(0,0), plotFrame(T,'g',.1);

odometry_state = zeros(6,0);
odometry_cov = ones(6,0);
pc=[];

state_mean = zeros(6,1);
state_covariance = eye(6)*1;
proc_noise = diag([.001 .001 .001 0.0001 0.0001 0.0001]);
obs_noise=diag([.01 .01 .01 .01 .01 .01]);

while(1)
	if index >last
		break;
	end
	
	% capturing frame
	[W,I]=mxNiPointCloud(context);
	dt=toc;tic;
	
	I = rgb2gray(I);
	% matching
	matcherMex('push',I);
	matcherMex('match',0);
	matches = matcherMex('getmatches',0);
	if size(matches,2)<15
		continue
	end
	%
	
	w0 = sparseKinect3D(matches(1:2,:),W0);
	w1 = sparseKinect3D(matches(3:4,:),W); 
	inl = logical( ( w0(3,:)>0 ) .* ( w1(3,:)>0 ));
	w0=w0(:,inl);
	w1=w1(:,inl);
	
	W0=W;
		
	% computing rigid transform
	Tr=rigidTransformation(w0,w1,'ransac');
	
	%storing the observation vector
	[roll,pitch,yaw]=dcm2angle(Tr(1:3,1:3));
	shift = Tr(1:3,4);
	ang = [roll pitch yaw]';
	obs_vector = [shift ; ang ];
	
	%ukf filtering
	[state_mean,state_covariance] ...
		= ukf_estimator1(...
					state_mean,...
					state_covariance,...
					eye(6),...
					proc_noise,...
					obs_vector,...
					eye(6),...
					obs_noise...
					);
			
						
	odometry_state(:,index)=state_mean;
	odometry_covariance(:,index)=diag(state_covariance);
	
	roll=odometry_state(4,index);
	pitch=odometry_state(5,index);
	yaw=odometry_state(6,index);
	Tr=[angle2dcm(yaw,pitch,roll) odometry_state(1:3,index);0 0 0 1];
	
	T=T*Tr;
	
	
	% plot statistics
	fprintf('matches: %i inliers: %.2f%s\n',size(matches,2),100*nnz(inl)/size(matches,2),char(37))
	subplot(3,3,[4 5 7 8]), scatter3(T(1,4),T(2,4),T(3,4),'r','+');%plotFrame(T,'r',.1)
	subplot(3,3,1), imshow(I);
	subplot(3,3,2), imshow(100*uint8(W(:,:,3)));
	subplot(3,3,3), plot(odometry_state(1:3,:)');
	subplot(3,3,6), plot(odometry_state(4:6,:)');

	drawnow
	
	
	index = index + step;% index increment

end
%%
%% Stop libvisio
matcherMex('close');

%% Stop the Kinect Process
mxNiDeleteContext(context);
clear context xml_config_file;

%%