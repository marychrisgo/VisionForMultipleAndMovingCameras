% Exercises based on the fundamental matrix
% Ex3: Computation of the essential matrix and euclidean reconstruction
% Uses ACT and extra_funs

clear, close all,

% include ACT_lite path
ACT_path = '../ACT_lite/';
addpath(genpath(ACT_path));
% include extra funs
extra_funs_path = '../extra_funs/';
addpath(genpath(extra_funs_path));

warning off
disp('************************************* START')

% load synthetic scene
load('data_F3lab.mat')

% Visualize the original scene (Euclidean)
figure, draw_scene(Q, K, R, t);
title('Original Scene. Ground truth')

% ------------------------------------------------------------------------
% 0. Select projected points (ideal or noisy version)
% ------------------------------------------------------------------------
% q_data = q;  % Uncomment this one to use ideal points
q_data = q_r; % Uncomment this one to use noisy points

% ------------------------------------------------------------------------
% 1. Compute the fundamental matrix and projection matrices. Use the
% MatFunProjectiveCalib function
% ------------------------------------------------------------------------
[F,P_est,Q_est,q_est] = MatFunProjectiveCalib(q_data);

disp(['Residual reprojection error. 8 point algorithm   = ' num2str( ErrorRetroproy(q_data,P_est,Q_est)/2 )]);
draw_reproj_error(q_data,P_est,Q_est);

% ------------------------------------------------------------------------
% 2. Visualization of ideal projectios, noisy projections and reprojected points
% ------------------------------------------------------------------------
for k=1:ncam
     figure()  
     hold on
     scatter(q(1,:,k),q(2,:,k),30,[1,0,0]);            % Ideal projections. Red
     scatter(q_r(1,:,k),q_r(2,:,k),30,[0,1,0]);        % Noisy projections. Green
     q_est = un_homogenize_coords(q_est);
     scatter(q_est(1,:,k),q_est(2,:,k),30,[0,0,1]);    % Reprojected points after projective calibration. Blue
     hold off
     title(sprintf('Image %d', k));
     % axis([-1000, 1000, -1000, 1000]);
     daspect([1, 1, 1]);
     pbaspect([1, 1, 1]);
end
hold off

% ------------------------------------------------------------------------
% 3. Obtain the essential matrix (E) from the fundamental matrix (F) and the
% intrinsic parameter matrices (K).
% ------------------------------------------------------------------------
E = K(:,:,2)' * F * K(:,:,1);

% ------------------------------------------------------------------------
% 4. Factorize the essential matrix with the 2 possible solutions for
% R. Use the function factorize_E to obtain R_est(:,:,1) and R_est(:,:,2) and T_est.
% ------------------------------------------------------------------------
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

% ------------------------------------------------------------------------
% 5. For each solution we obtain an Euclidean solution and we visualize it.
% ------------------------------------------------------------------------
npoints = size(q_data,2);
Q_euc = zeros(4,npoints,2); % Variable for recontructed points
P_euc = zeros(3,4,2);       % Variable for projection matrices
figNo=figure;

for sol=1:4
    % Euclidean triangulation to obtain the 3D points (use TriangEuc)
    Q_euc = TriangEuc(Rcam(:,:,2,sol),Tcam(:,2,sol),K,q_est);
       
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

disp('************************************* END')

