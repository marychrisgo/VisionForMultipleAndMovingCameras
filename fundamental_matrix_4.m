% Exercises based on the fundamental matrix
% Ex4: Computation of fundamental matrix and affine reconstruction from
% vanishing points (points at infinity). Eucludean reconstruction from the
% knowledge of the scene structure

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
load('data_F4lab.mat')

draw_scene(Q, K, R, t);
draw_3D_cube_segments(Q);

% ------------------------------------------------------------------------
% 0. Select projected points (ideal or noisy version)
% ------------------------------------------------------------------------
% q_data = q;  % Uncomment this one to use ideal points
q_data = q_r; % Uncomment this one to use noisy points

% ------------------------------------------------------------------------
% 1. Compute the fundamental matrix and projection matrices. Use the
% MatFunProjectiveCalib function
% ------------------------------------------------------------------------
[F_est,P_est,Q_est,q_est] = MatFunProjectiveCalib(q_data);

disp(['Resudual reprojection error. 8 point algorithm   = ' num2str( ErrorRetroproy(q_data,P_est,Q_est)/2 )]);
draw_reproj_error(q_data,P_est,Q_est);

% ------------------------------------------------------------------------
% 2. Visualize results. Projected points and projective reconstruction
% ------------------------------------------------------------------------
% Projected points
draw_projected_cube(q_est)

% Projective reconstruction. Use draw_3Dcube 
draw_3Dcube(Q_est);

% ------------------------------------------------------------------------
% 3. Compute the affine reconstrucion based on vanishing points (projections of points at infinity on your images)
% ------------------------------------------------------------------------
% Compute the projection of the 3 vanishing points in each image. Use cross() (cross-vector product)
% The 8 first points in q_data are the positions of the cube vertices. 
% Calculate crossing points of parallel lines (in space) from the vertices of that cube 
% ...

% It is useful to check that the computed vanishing points are correct. You
% can use the draw_lines() for that purpose
for k=1:ncam
    figure();
    hold on
    scatter(q(1,:,k),q(2,:,k),30,[1,0,0]);
    draw_proj_cube(q(:,:,k));
    % draw the projection of lines and vanishing points in each image
    % ...
    % ...

end

% Build a (3,3,ncam) matrix with the position (hom. coords) of the vanishing points in each image
% ...

% Obtain the 3D points at the infinity plane by linear triangulation. Use
% the linear_triangulation() function
% ...

% The following function returns the hom coords. of a plane if X(3,4) contains 3 points of that plane 
% plane = NumKernel(X.'); ...

% Build the P3 transformation H_aff that transforms the projective reconstruction into a affine reconstruction
% H_aff = ....

% Apply H_aff to points and projection cameras
Qa = H_aff * Q_est;
Pa = zeros(3, 4, 2);
Pa(:, :, 1) = P_est(:, :, 1) * inv(H_aff);
Pa(:, :, 2) = P_est(:, :, 2) * inv(H_aff);
Xa = H_aff * X;


% ------------------------------------------------------------------------
% 4. % Visualize affine reconstruction. Use draw_3Dcube() 
% ------------------------------------------------------------------------
draw_3Dcube(Qa);

% ------------------------------------------------------------------------
% 5. Compute the Euclidean reconstruction, knowing that we have a cubic
% scene
% ------------------------------------------------------------------------
% Use calc_reference_homography(Q): Computes the space homography Heuc that transforms the canonical reference 
% (1 0 0 0) (0 1 0 0) (0 0 1 0) (0 0 0 1) (1 1 1 1) into the 5 3D points that are stored in Q 

Heuc = inv(calc_reference_homography(Q_ref));

% Apply that transformation to the projective or affine points and
% matrices. Depending on which you used to obtain the transformation
Q_euc = Heuc * Qa;
P_euc = zeros(3, 4, 2);
P_euc(:, :, 1) = Pa(:, :, 1) * inv(Heuc);
P_euc(:, :, 2) = Pa(:, :, 2) * inv(Heuc);


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

disp('************************************* END')





