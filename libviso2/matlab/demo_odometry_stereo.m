% demonstrates stereo visual odometry on synthetic test data
disp('===========================');
clear all; close all; dbstop error;

% bucketing parameters
max_features         = int32(200);
bucket_width         = 50;
bucket_height        = 50;

% odometry/camera parameters
f                    = 340;
cu                   = 330;
cv                   = 94;
base                 = 0.6;
varMeasurements      = 0.2;
deltaT               = 0.1;

intrinsicCalibration = [f 0 cu;0 f cv;0 0 1];
extrinsicRotation    = eye(3);
extrinsicTranslation = [-base;0;0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     create test trajectory      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=0;

% create random stream with seed
rstream1=RandStream('mrg32k3a','Seed',1+k); % trajectory
rstream2=RandStream('mrg32k3a','Seed',2+k); % features

% synthesize transformation matrices for all frames
alpha = 0;
for i = 1:100
  alpha = 0.95*(alpha + 0.03*randn(rstream1));
  Tr_ref{i} = [cos(alpha) 0 sin(alpha) 0;0 1 0 0; -sin(alpha) 0 cos(alpha) -2; 0 0 0 1];
  if i==1
    Tr_ref_total{i} = Tr_ref{i}^-1;
  else
    Tr_ref_total{i} = Tr_ref_total{i-1}*Tr_ref{i}^-1;
  end
end

% projection matrix from left to right camera coordinate system
Tr_left_right = [extrinsicRotation,extrinsicTranslation;0 0 0 1];

% create 3d points and 2d matches for all frames
num_features = 50;
for i=1:length(Tr_ref)
  
  % create 3d reference points
  p3_ref = [];
  for j=1:num_features
    p3_ref = [p3_ref; randn(rstream2)*10 randn(rstream2)*5 15+rand(rstream2)*20];
  end
  
  % 3d points in camera coordinate systems
  p3{i}.lp = p3_ref;
  p3{i}.rp = project(p3_ref,Tr_left_right);
  p3{i}.lc = project(p3_ref,Tr_ref{i});
  p3{i}.rc = project(p3_ref,Tr_left_right*Tr_ref{i});
  
  % project to image planes
  p2{i}.lp = project(p3{i}.lp,intrinsicCalibration)+randn(rstream2,size(p3{i}.lp,1),2)*sqrt(varMeasurements);
  p2{i}.rp = project(p3{i}.rp,intrinsicCalibration)+randn(rstream2,size(p3{i}.rp,1),2)*sqrt(varMeasurements);
  p2{i}.lc = project(p3{i}.lc,intrinsicCalibration)+randn(rstream2,size(p3{i}.lc,1),2)*sqrt(varMeasurements);
  p2{i}.rc = project(p3{i}.rc,intrinsicCalibration)+randn(rstream2,size(p3{i}.rc,1),2)*sqrt(varMeasurements);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     compute visual odometry     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% init matcher + odometry object
visualOdometryMex('init',f,cu,cv,base);
                       
% create figure
figure; axis equal; hold on;
plot(0,0,'-xr','LineWidth',1);
plot(0,0,'-xb','LineWidth',1);
legend('Ground truth','Visual odometry');

% for all frames do
for i=1:length(Tr_ref)

  % estimate egomotion
  p_matched = [p2{i}.lp p2{i}.rp p2{i}.lc p2{i}.rc]';
  success = visualOdometryMex('update',deltaT,p_matched);

  % get state vector and inliers
  Tr_odo{i} = visualOdometryMex('gettransformation');
  inliers   = visualOdometryMex('getinliers');
  
  % accumulate total motion
  if i==1
    Tr_odo_total{i} = Tr_odo{i}^-1;
  else
    Tr_odo_total{i} = Tr_odo_total{i-1}*Tr_odo{i}^-1;
    
    % show resulting trajectory
    plot([Tr_ref_total{i-1}(1,4) Tr_ref_total{i}(1,4)], ...
         [Tr_ref_total{i-1}(3,4) Tr_ref_total{i}(3,4)],'-xr','LineWidth',1);
    plot([Tr_odo_total{i-1}(1,4) Tr_odo_total{i}(1,4)], ...
         [Tr_odo_total{i-1}(3,4) Tr_odo_total{i}(3,4)],'-xb','LineWidth',1);
       
    % redraw
    pause(0.01);
    refresh;
  end
  
  % output
  disp(['Frame: ' num2str(i) ', Inliers: ' num2str(length(inliers))]);
end

% release objects
visualOdometryMex('close');
