% takes a monocular input sequence and computes 3d point clouds
disp('===========================');
clear all; close all;

% image data
% (download sequence '2010_03_09_drive_0019' from www.cvlibs.net)
img_dir     = 'C:\Users\Redjan\Desktop\workspace\data\outdoor\2010_03_09_drive_0019';
first_frame = 0;
last_frame  = 300;

% matching parameters
nms_n                  = 3;   % non-max-suppression: min. distance between maxima (in pixels)
nms_tau                = 50;  % non-max-suppression: interest point peakiness threshold
match_binsize          = 50;  % matching bin width/height (affects efficiency only)
match_radius           = 200; % matching radius (du/dv in pixels)
match_disp_tolerance   = 1;   % du tolerance for stereo matches (in pixels)
outlier_disp_tolerance = 5;   % outlier removal: disparity tolerance (in pixels)
outlier_flow_tolerance = 5;   % outlier removal: flow tolerance (in pixels)
multi_stage            = 1;   % 0=disabled,1=multistage matching (denser and faster)
half_resolution        = 0;   % 0=disabled,1=match at half resolution, refine at full resolution
refinement             = 2;   % refinement (0=none,1=pixel,2=subpixel)

% bucketing parameters
max_features         = 2;
bucket_width         = 50;
bucket_height        = 50;
deltaT               = 0.1;

% init matcher + odometry objects
matcherMex('init',nms_n,nms_tau,match_binsize,match_radius,match_disp_tolerance,outlier_disp_tolerance,outlier_flow_tolerance,multi_stage,half_resolution,refinement);
visualOdometryMex('init',679,660,187,1.6,-0.08);
reconstructionMex('init',679,660,187);

% init transformation matrix array
Tr_total{1} = eye(4);

% create figure
figure('Color',[1 1 1]),axes('Position',[0,0,1,1]); axis equal, axis off;

% for all frames do
failure = 0;
for frame=first_frame:last_frame
  
  % read current images
  I1 = keyFinder(frame);

  % keep previous image, if visual odometry estimate failed
  if failure==0
    matcherMex('push',I1');
  else
    matcherMex('replace',I1');
  end

  % start matching after reading 2nd frame  
  if frame>1

    % match features and also get a bucketed version
    matcherMex('match',0);
    p_matched_full = matcherMex('getmatches',0);
    i_matched_full = matcherMex('getidx',0);
    matcherMex('bucketing',max_features,bucket_width,bucket_height);
    p_matched = matcherMex('getmatches',0);

    % estimate egomotion, failure+1 is the number of frames between matches
    success = visualOdometryMex('update',(failure+1)*deltaT,p_matched);
    inliers = visualOdometryMex('getinliers');
    
    % if visual odometry estimate was successful
    if success==1
      
      % get transformation matrix
      Tr = visualOdometryMex('gettransformation');
      
      % reconstruct non-road points of finished tracks
      reconstructionMex('update',p_matched_full,i_matched_full,Tr,2,3,30,3);
      failure = 0;
      
    % otherwise keep current position
    % (this happens if motion is very small)
    else
      Tr      = eye(4);
      failure = failure+1;
    end
    
    % accumulate total transformation
    Tr_total{frame} = Tr_total{frame-1}*Tr^-1;
    
    % get reconstructed 3d points
    p = reconstructionMex('getpoints');
    
    % update point cloud and draw trajectory
    cla; hold on; view(0,0);
    cla,plot3(p(1,:),p(2,:),p(3,:),'.b','MarkerSize',1),axis equal;
    for i=2:frame
      plot3([Tr_total{i-1}(1,4) Tr_total{i}(1,4)], ...
            [Tr_total{i-1}(2,4) Tr_total{i}(2,4)], ...
            [Tr_total{i-1}(3,4) Tr_total{i}(3,4)],'-xr','LineWidth',1);
    end
    pause(0.05);
    refresh;

    % status
    disp(['Frame: ' num2str(frame), ' Matches: ' num2str(size(p_matched,2)), ' Inliers: ' num2str(100*length(inliers)/size(p_matched,2),'%3.1f') '%']);
  end
end

% release objects
matcherMex('close');
visualOdometryMex('close');
reconstructionMex('close');
