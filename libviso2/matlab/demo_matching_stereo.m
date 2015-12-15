% demonstrates sparse stereo
disp('===========================');
clear all; dbstop error; close all;

% read images from file
I1c = imread('../img/I1c.png');
I2c = imread('../img/I2c.png');

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
refinement             = 1;   % refinement (0=none,1=pixel,2=subpixel)

% init matcher
matcherMex('init',nms_n,nms_tau,match_binsize,match_radius,match_disp_tolerance,outlier_disp_tolerance,outlier_flow_tolerance,multi_stage,half_resolution,refinement);

% push back images
tic
matcherMex('push',I1c',I2c');
disp(['Feature detection: ' num2str(toc) ' seconds']);

% match images
tic
matcherMex('match',1);
p_matched = matcherMex('getmatches',1);
disp(['Feature matching:  ' num2str(toc) ' seconds']);

% close matcher
matcherMex('close');

% show matching results
disp(['Number of matched points: ' num2str(length(p_matched))]);
disp('Plotting ...');
plotMatch(I1c,p_matched,1);
