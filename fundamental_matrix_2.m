% Exercises based on the fundamental matrix
% Ex2: Generate a synthetic scene and 2 cameras. Then compute F from image correspondences
% Obtain P up to a projective transformantion and check that the
% reconstructed 3D points points are different from the original scene, but projections (images), are preserved
% Uses the ACT toolbox

clear all, close all,

% include ACT_lite path
ACT_path = '../ACT_lite/';
addpath(genpath(ACT_path));
% include extra funs
extra_funs_path = '../extra_funs/';
addpath(genpath(extra_funs_path));

warning off
disp('************************************* START')

% load synthetic scene
load('data_F2lab.mat')

% Visualize the original scene (Euclidean)
figure, draw_scene(Q, K, R, t);
title('Original Scene. Ground truth')

% visualize projected points. Original and noisy version
draw_projected_points_noise(q,q_r)

% ------------------------------------------------------------------------
% 0. Select projected points (ideal or noisy version)
% ------------------------------------------------------------------------
q_data = q;  % Uncomment this one to use ideal points
% q_data = q_r; % Uncomment this one to use noisy points

% ------------------------------------------------------------------------
% 1. Compute the fundamental matrix and projection matrices. Use MatFunProjectiveCalib function 
% ------------------------------------------------------------------------
[F,P_est,Q_est,q_est] = MatFunProjectiveCalib(q_data);

q_est = un_homogenize_coords(q_est);
disp(['Residual reprojection error. 8 point algorithm   = ' num2str( ErrorRetroproy(q_data,P_est,Q_est)/2 )]);
draw_reproj_error(q_data,P_est,Q_est);

% ------------------------------------------------------------------------
% 2. Visualize results
% ------------------------------------------------------------------------
% Compare data points and reprojected points
draw_projected_points_noise(q_data,q_est)

% ------------------------------------------------------------------------
% 3. Compute and visualize epipolar lines
% ------------------------------------------------------------------------
% ep_lines(3,npoints,ncam)
ep_lines = zeros(3,npoints,ncam);

for i = 1:npoints
    ep_lines(:,i,1) = F' * q_est(:,i,2);
    ep_lines(:,i,2) = F * q_est(:,i,1);
end

for k=1:ncam
    xmin = min(q_est(1,:,k));
    xmax = max(q_est(1,:,k));
    ymin = min(q_est(2,:,k));
    ymax = max(q_est(2,:,k));
    figure();
    hold on
    scatter(q_est(1,:,k),q_est(2,:,k),30,[1,0,0]);
    draw_lines(ep_lines(:,:,k),xmin,xmax,ymin,ymax,'b');
end

% Projective reconstruction. Use draw_3Dcube 
draw_3Dcube(Q_est);

