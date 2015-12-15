% demonstrates monocular feature tracking (via feature indices)
clear all; dbstop error; close all;

% matching parameters
nms_n                  = 3;   % non-max-suppression: min. distance between maxima (in pixels)
nms_tau                = 50;  % non-max-suppression: interest point peakiness threshold
match_binsize          = 50;  % matching bin width/height (affects efficiency only)
match_radius           = 100; % matching radius (du/dv in pixels)
match_disp_tolerance   = 1;   % du tolerance for stereo matches (in pixels)
outlier_disp_tolerance = 5;   % outlier removal: disparity tolerance (in pixels)
outlier_flow_tolerance = 5;   % outlier removal: flow tolerance (in pixels)
multi_stage            = 1;   % 0=disabled,1=multistage matching (denser and faster)
half_resolution        = 0;   % 0=disabled,1=match at half resolution, refine at full resolution
refinement             = 1;   % refinement (0=none,1=pixel,2=subpixel)

% init matcher
matcherMex('init',nms_n,nms_tau,match_binsize,match_radius,match_disp_tolerance,outlier_disp_tolerance,outlier_flow_tolerance,multi_stage,half_resolution,refinement);

% push back first image
I = imread('../img/I1_000000.png');
matcherMex('push',I');

% feature tracks
tracks = {};

% start matching
for i=1:6
  I = imread(['../img/I1_' num2str(i,'%06d') '.png']);
  tic
  matcherMex('push',I');
  disp(['Feature detection: ' num2str(toc) ' seconds']);
  tic
  matcherMex('match',0);
  p_matched{i} = matcherMex('getmatches',0);
  i_matched{i} = matcherMex('getidx',0);
  disp(['Feature matching:  ' num2str(toc) ' seconds']);
end

% close matcher
matcherMex('close');

% show matching results
disp('Plotting ...');
plotTrack(I,p_matched,i_matched);
