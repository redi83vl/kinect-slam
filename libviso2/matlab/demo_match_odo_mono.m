% demonstrates monocular visual odometry on a real sequence
% note: a lot of time is spent on drawing feature matches
disp('===========================');
clear all; close all; dbstop error;

% image data
% (download sequence '2010_03_09_drive_0019' from www.cvlibs.net)
img_dir     = '/home/geiger/5_Data/karlsruhe_dataset/2010_03_09_drive_0019';
first_frame = 0;
last_frame  = 372;

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
refinement             = 1;   % refinement (0=none,1=pixel,2=subpixel)

% bucketing parameters
max_features         = 10;
bucket_width         = 50;
bucket_height        = 50;
deltaT               = 0.1;

% init matcher + odometry objects
matcherMex('init',nms_n,nms_tau,match_binsize,match_radius,match_disp_tolerance,outlier_disp_tolerance,outlier_flow_tolerance,multi_stage,half_resolution,refinement);
visualOdometryMex('init',679,660,187,1.6,-0.08);

% init transformation matrix array
Tr_total{1} = eye(4);

% create figure
figure('Color',[1 1 1]);
subplot(2,1,2); axis equal, grid on;
set(gca,'XTick',-500:10:500);
set(gca,'YTick',-500:10:500);
set(gca,'ZTick',-500:10:500);

% for all frames do
failure = 0;
k = 2;
for frame=first_frame:last_frame
  
  % read current images
  I1 = imread([img_dir '/I1_' num2str(frame,'%06d') '.png']);

  % keep previous image, if visual odometry estimate failed
  if failure==0
    matcherMex('push',I1');
  else
    matcherMex('replace',I1');
  end

  % start matching after reading 2nd frame  
  if frame>first_frame+1

    % match features and also get a bucketed version
    matcherMex('match',0);
    p_matched_full = matcherMex('getmatches',0);
    matcherMex('bucketing',max_features,bucket_width,bucket_height);
    p_matched = matcherMex('getmatches',0);

    % estimate egomotion, failure+1 is the number of frames between matches
    success = visualOdometryMex('update',(failure+1)*deltaT,p_matched);
    inliers = visualOdometryMex('getinliers');
    
    % if visual odometry estimate was successful
    if success==1
      
      % get transformation matrix
      Tr = visualOdometryMex('gettransformation');
      failure = 0;
      
    % otherwise keep current position
    % (this happens if motion is very small)
    else
      Tr      = eye(4);
      failure = failure+1;
    end
    
    % accumulate total transformation
    Tr_total{k} = Tr_total{k-1}*Tr^-1;
    k=k+1;
    
    % plot match
    subplot(2,1,1);
    plotMatch(I1,p_matched,0,inliers);
    
    % plot trajectory
    subplot(2,1,2);
    cla; hold on; view(0,0);
    for i=1:length(Tr_total)-1
      plot3([Tr_total{i}(1,4) Tr_total{i+1}(1,4)], ...
            [Tr_total{i}(2,4) Tr_total{i+1}(2,4)], ...
            [Tr_total{i}(3,4) Tr_total{i+1}(3,4)],'-xr','LineWidth',1);
    end
    
    % update view
    pause(0.1);
    refresh;

    % status
    disp(['Frame: ' num2str(frame), ' Matches: ' num2str(size(p_matched,2)), ' Inliers: ' num2str(100*length(inliers)/size(p_matched,2),'%3.1f') '%']);
  end
end

% release objects
matcherMex('close');
visualOdometryMex('close');
