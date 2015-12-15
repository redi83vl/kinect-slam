% demonstrates monocular visual odometry on synthetic test data
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
varMeasurements      = 0.2;
cam_height           = 1.6;
deltaT               = 0.1;

intrinsicCalibration = [f 0 cu;0 f cv;0 0 1];

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

% create 3d points and 2d matches for all frames
num_features = 100;
for i=1:length(Tr_ref)
  
  p3_ref = [];
  
  % create random 3d reference points
  for j=1:num_features
    p3_ref = [p3_ref; randn(rstream2)*10 randn(rstream2)*5 10+rand(rstream2)*20];
  end
  
  % create 3d reference points on ground plane
  for j=1:num_features
    p3_ref = [p3_ref; randn(rstream2)*10 cam_height 10+rand(rstream2)*20];
  end
  
  % 3d points in camera coordinate systems
  p3{i}.lp = p3_ref;
  p3{i}.lc = project(p3_ref,Tr_ref{i});
  
  % project to image planes
  p2{i}.lp = project(p3{i}.lp,intrinsicCalibration)+randn(rstream2,size(p3{i}.lp,1),2)*sqrt(varMeasurements);
  p2{i}.lc = project(p3{i}.lc,intrinsicCalibration)+randn(rstream2,size(p3{i}.lc,1),2)*sqrt(varMeasurements);
  
  % outliers
  for j=1:5:size(p2{i}.lp,1)
    p2{i}.lp(j,:) = p2{i}.lp(j,:)+randn(rstream2)*5;
    p2{i}.lc(j,:) = p2{i}.lc(j,:)+randn(rstream2)*5;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     compute visual odometry     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% init matcher + odometry object
visualOdometryMex('init',f,cu,cv,cam_height);

% create figure
figure; axis equal; hold on;
plot(0,0,'-xr','LineWidth',1);
plot(0,0,'-xb','LineWidth',1);
legend('Ground truth','Visual odometry');

% for all frames do
for i=1:length(Tr_ref)
  
  % estimate egomotion
  p_matched = [p2{i}.lp p2{i}.lc]';
  success = visualOdometryMex('update',deltaT,p_matched);

  % get state vector and inliers
  Tr_odo{i} = visualOdometryMex('gettransformation');
  
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
end

% release objects
visualOdometryMex('close');
