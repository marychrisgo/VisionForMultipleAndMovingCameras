clear all
close all
clc

% include ACT_lite path
ACT_path = './ACT_lite/';
addpath(genpath(ACT_path));
% include extra funs
extra_funs_path = './extra_funs/';
addpath(genpath(extra_funs_path));

%% PARAMETER DEFINITION %%
%%define dataset path
params.Directory    = fullfile('Dataset4/');

%%Detector parameters
params.detector     =  'SIFT'; %'LoG_ss', 'SURF', 'SIFT','DoG_ss', 'KAZE' 'DoH' ...
params.nscales      =       15;
params.noctaves     =        8;
params.sigma0       =      1.6; % as we are using Matlab functions this is the minimum value allowed
params.npoints      =      500;

%%Descriptor parameters (see doc extractFeatures for explanation and additional parameters)
params.descriptor   =  'DSP-SIFT'; % 'SIFT', 'SURF', 'DSP-SIFT', 'KAZE'...
params.desOnDecom   =   false; % describe on scale-space (linear or non-linear) decomposition (if available)
params.Upright      =   false; % set to true to avoid orientation estimation.
% for DSP-SIFT
params.dsp.ns       =      10;% number of sampled scales
params.dsp.sc_min   =     1/6;% smallest scale (relative to detection)
params.dsp.sc_max   =       3;% largest scale (relative to detection);    

%%Matching parameters (see doc matchFeatures for explanation and additional parameters)
params.MaxRatio     =   0.6;  %= 
params.Metric       =  'SSD';
%% END OF PARAMETER DEFINITION %%

%% addpaths
addpath(genpath('./detectors/'));
addpath(genpath('./descriptors/'));
addpath(genpath('./toolbox/'));

%% preload dataset
params.Scene = imageDatastore(params.Directory);
numImages    = numel(params.Scene.Files);

%% initialize (sort of)
ima{numImages}           = [];
points{numImages}        = [];
decomposition{numImages} = [];
features{numImages}      = [];

%% get sigmas
k = 1;
params.sigmas = zeros(1,params.noctaves*params.nscales);
for o = 0:params.noctaves-1
    params.sigmas(k:(k+params.nscales-1)) = params.sigma0.*pow2([0:(params.nscales-1)]/params.nscales + o);
    k = k+params.nscales;
end

%% detect & describe
for j = 1:numImages
%% Load and convert images %%
SCALE = 0.4
ima{j}      =    imresize(readimage(params.Scene, j), SCALE);
gima        =      im2double(rgb2gray(ima{j})); 

%% PoI Detection %%
[points{j},decomposition{j}] =  myDetector(gima,params);

%% PoI Description %%
[features{j},points{j}]      =  myDescriptor(points{j},decomposition{j},params);

end

%load('point.mat');
%load('features.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. compute consistent point matches from N views
point_matrix = n_view_matching(points, features, ima, params.MaxRatio, params.Metric);
npoints = size(point_matrix,2);

homography = ones(3,npoints,numImages);
homography(1:2,:,:) = point_matrix;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. compute Fundamental matrix and in initial projective reconstruction
camera(:,:,1)= homography(:,:,1); 
camera(:,:,2)= homography(:,:,3);

% mean re-projection error and the reprojection error histogram
[F, P_est2,Q_est2,q_est2] = MatFunProjectiveCalib(camera);

vgg_gui_F(ima{1},ima{numImages},F');
%disp(['Residual reprojection error (cameras)  = ' num2str( ErrorRetroproy(camera,P_est2,Q_est2)/2 )]);
%draw_reproj_error(camera,P_est2,Q_est2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Projection matrix of the rest of the cameras
n_camera = zeros(3,4,numImages);
n_camera(:,:,1) = P_est2(:,:,1);
n_camera(:,:,numImages) = P_est2(:,:,2);
for i=1:numImages
    n_camera(:,:,i) = PDLT_NA(homography(:,:,i),Q_est2,0,0); % 2D, 3D,SVD, does not print costs in the command line
end

% statistics
%disp(['Residual reprojection after resectioning  = ' num2str( ErrorRetroproy(homography,n_camera,Q_est2)/2 )]);
%draw_reproj_error(homography,n_camera,Q_est2);

vp = ones(npoints,numImages);
[projection_BA,Q_BA,reprojected_BA] = BAProjectiveCalib(homography,n_camera,Q_est2,vp);

% statistics
%disp(['Residual reprojection after Projective Bundle Adjustment = ' num2str( ErrorRetroproy(homography,projection_BA,Q_BA)/2)]);
%draw_reproj_error(homography,projection_BA,Q_BA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Re-compute F
recomputed_F = vgg_F_from_P(projection_BA(:,:,2),projection_BA(:,:,3));
vgg_gui_F(ima{1},ima{numImages},recomputed_F');

% 5. Euclidean reconstruction
reprojected_BA2(:,:,1) = reprojected_BA(:,:,1);
reprojected_BA2(:,:,2) = reprojected_BA(:,:,2);

load('A_matrix.mat');
A_small = A_small * SCALE
% A_small(1,2) = -A_small(1,2)* cos(A_small(1,2)); % or cos 90?
% A_small(2,2) = A_small(1,2)/sin(A_small(1,2));
A_small(1,2) = acos(-A_small(1,2)/A_small(1,1)); % or cos 90?
A_small(2,2) = A_small(2,2)/sin(A_small(1,2));

% intrinsic parameter
K(:,:,1) = A_small; 
K(:,:,2) = A_small;

E = K(:,:,2)' * recomputed_F * K(:,:,1);
[R_est,T_est] = factorize_E(E);

% ------------------------------------------------------------------------
% Save the 4 solutions (R,t) in the structures Rcam(3,3,cam,sol), T(3,cam,sol),
% where cam indicates the camera number and sol indicates the solution number (1, 2, 3 or 4).
% ------------------------------------------------------------------------
Rcam = zeros(3,3,2,4);
Tcam = zeros(3,2,4);

% first camera
Rcam(:,:,1,1) = eye(3,3);
Rcam(:,:,1,2) = eye(3,3);
Rcam(:,:,1,3) = eye(3,3);
Rcam(:,:,1,4) = eye(3,3);

% Tcam of 1 is null translation

% second camera
Rcam(:,:,2,1) = R_est(:,:,1);
Rcam(:,:,2,2) = R_est(:,:,1);
Rcam(:,:,2,3) = R_est(:,:,2);
Rcam(:,:,2,4) = R_est(:,:,2);

Tcam(:,2,1) = T_est;
Tcam(:,2,2) = -T_est;
Tcam(:,2,3) = T_est;
Tcam(:,2,4) = -T_est;

npoints = size(camera,2);
Q_euc = zeros(4,npoints,4); % Variable for recontructed points
P_euc = zeros(3,4,2,4);     % Variable for projection matrices
figNo=figure;
reprojectedPoints = zeros(3,npoints,2,4);

for sol=1:4
    % Euclidean triangulation to obtain the 3D points (use TriangEuc)
    Q_euc = TriangEuc(Rcam(:,:,2,sol),Tcam(:,2,sol),K,reprojected_BA2);
    Q_euc2(:,:,sol) = Q_euc;
    % visualize 3D reconstruction
    figure();
    draw_scene(Q_euc, K, Rcam(:,:,:,sol), Tcam(:,:,sol));
    title(sprintf('Solution %d', sol));
     
    % Compute the projection matrices from K, Rcam, Tcam
    for k=1:2
       P_euc(:,:,k) = K(:,:,k) * [Rcam(:,:,k,sol), -Rcam(:,:,k,sol) * Tcam(:,k,sol)]; %euclidean
    end
    
    % Obtain the re-projected points q_rep
    q_rep = zeros(3, npoints, 2);
    q_rep(:,:,1) = P_euc(:,:,1)* Q_euc;
    q_rep(:,:,2) = P_euc(:,:,2)* Q_euc;
    reprojectedPoints(:,:,:,sol) = q_rep;
    
    % Visualize reprojectd points to check that all solutions correspond to
    % the projected images
    q_rep = un_homogenize_coords(q_rep);
    for k=1:2
      figure(figNo); subplot(4,2,2*(sol-1)+k); scatter(q_rep(1,:,k),q_rep(2,:,k),30,[1,0,0]);
      title(sprintf('Reprojection %d, image %d', sol, k));
      daspect([1, 1, 1]);
      pbaspect([1, 1, 1]);
      axis([-1000, 1000, -1000, 1000]);
    end
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------------------------------------------------------------------------
% 6. Compute the intrinsic and extrinsic parameters from the euclidean
% projection matrices. Use CameraMatrix2KRC
% ------------------------------------------------------------------------
[K_euc,R_euc,C_euc] = CameraMatrix2KRC(P_euc);

% Visualize results. Check both possible solutions
draw_3Dcube(Q_euc);
figure;
draw_scene(Q_euc, K_euc, R_euc, C_euc(1:3,:));
draw_3D_cube_segments(Q_euc);
figure;
S=-eye(4); S(4,4)=1;
draw_scene(S*Q_euc, K_euc, R_euc, -C_euc(1:3,:));
draw_3D_cube_segments(S*Q_euc);


disp(['Residual reprojection after Euclidean Reconstruction = ' num2str( ErrorRetroproy(homography,n_camera,Q_est2)/2)]);
draw_reproj_error(homography,n_camera,Q_est2);

disp('************************************* END')
