% Exercises for N views
% N view reconstruction. Initial reconstruction using the F matrix
% between 2 cameras. Resection for the rest of cameras and bundle adjustment
% Euclidean reconstruction is obtained using the knowledge of
% the inf plane and scene structure

clear, close all,

% include ACT_lite path
ACT_path = './ACT_lite/';
addpath(genpath(ACT_path));
% include extra funs
extra_funs_path = './extra_funs/';
addpath(genpath(extra_funs_path));

warning off
disp('************************************* START')

load('data_BA_lab.mat')

% ------------------------------------------------------------------------
% 0. Select projected points (noisy version)
% ------------------------------------------------------------------------
q_data = q_r; % Uncomment this one to use noisy points

% ------------------------------------------------------------------------
% 1. Visualize results. Compare ideal (q) and noisy projections (q_r).
% ------------------------------------------------------------------------
% Projected points
draw_scene(Q, K(:,:,:), R(:,:,:), t(:,:));
draw_3D_cube_segments(Q);
title('Ground truth scene')

draw_projected_cube_noise(q,q_data); 


% ------------------------------------------------------------------------
% 2. Compute the fundamental matrix using the first and last cameras
% of the camera set (N cameras) # selecting 2 cameras out of
% 9
% ------------------------------------------------------------------------
q_2cams(:,:,1)=q_data(:,:,1); 
q_2cams(:,:,2)=q_data(:,:,ncam);

[F, P_2cam_est,Q_2cam_est,q_2cam_est] = MatFunProjectiveCalib(q_2cams);

disp(['Mean re-projection error of initial reconstruction (2 cams)   = ' num2str( ErrorRetroproy(q_2cams,P_2cam_est,Q_2cam_est)/2 )]);
draw_reproj_error(q_2cams,P_2cam_est,Q_2cam_est);

% ------------------------------------------------------------------------
% 3. Resection. Obtain the projection matrices of the rest of cameras using the PDLT_NA function 
% ------------------------------------------------------------------------

P_rep = zeros(3,4,ncam); % projection matrix is 3x4
P_rep(:,:,[1,ncam]) = P_2cam_est;
   
for i = 2:ncam-1
    P_rep(:,:,i) =  PDLT_NA(q_data(:,:,i), Q_2cam_est);
end

% ------------------------------------------------------------------------
% 4. Compute the statistics of the reprojection error for the initial projective reconstruction
% ------------------------------------------------------------------------
disp(['Mean re-projection error of initial reconstruction (9 cams)   = ' num2str( ErrorRetroproy(q_data,P_rep,Q_2cam_est)/2 )]);
draw_reproj_error(q_data,P_rep,Q_2cam_est);

% ------------------------------------------------------------------------
% 5. Proceed as in previous lab to extract the metric reconstruction from
% the initial N-view reconstruction
% ------------------------------------------------------------------------

p_l = zeros(3,6,ncam);
    % we need three pairs of parallel lines
for cam = 1:ncam
    p_l(:,1,cam) =cross(q_data(:,1,cam), q_data(:,5,cam));
    p_l(:,2,cam) =cross(q_data(:,2,cam), q_data(:,6,cam));
    p_l(:,3,cam) =cross(q_data(:,1,cam), q_data(:,2,cam));
    p_l(:,4,cam) =cross(q_data(:,5,cam), q_data(:,6,cam));
    p_l(:,5,cam) =cross(q_data(:,1,cam), q_data(:,4,cam));
    p_l(:,6,cam) =cross(q_data(:,5,cam), q_data(:,8,cam));
end

v_p = zeros(3,3,ncam); 
for cam =1:ncam
    v_p(:,1,cam) = cross(p_l(:,1,cam),p_l(:,2,cam));
    v_p(:,2,cam) = cross(p_l(:,3,cam),p_l(:,4,cam));
    v_p(:,3,cam) = cross(p_l(:,5,cam),p_l(:,6,cam));
    v_p(:,1,cam) =v_p(:,1,cam)/v_p(3,1,cam);
    v_p(:,2,cam) =v_p(:,2,cam)/v_p(3,2,cam);
    v_p(:,3,cam) =v_p(:,3,cam)/v_p(3,3,cam);
end
for k=1:ncam
    figure();
    hold on
    scatter(q(1,:,k),q(2,:,k),30,[1,0,0]);
    draw_proj_cube(q(:,:,k));
    scatter(v_p(1,:,k),v_p(2,:,k),'r');
    draw_lines(p_l(:,:,k),-8000,10000,-8000,8000,'b')
    % ...

end

v_p_3d = linear_triangulation(v_p,P_rep);
plane = NumKernel(v_p_3d.'); ...
H_aff = [eye(3,3) zeros(3,1); plane'];

Qa = H_aff*Q_2cam_est;
Pa(:,:,1) = P_rep(:,:,1)*inv(H_aff);
Pa(:,:,2) = P_rep(:,:,2)*inv(H_aff);
draw_3Dcube(Qa);

%Euclidean reconstruction

Q5 = zeros(4,5); % we need 5 points, the first three should be vanishing points, then two opposite vertices for the cube
Q5(:,1:3) = v_p_3d; %4*3 matrix
Q5(:,4)= Q_2cam_est(:,1); %4*1 coords for the 1st vertice 
Q5(:,5)= Q_2cam_est(:,7); %4*1 coords for the 7th vertice
Heuc = calc_reference_homography(Q5);

Q_euc = inv(Heuc)*Q_2cam_est;
for cam = 1:ncam
    P_euc(:,:,cam) = P_rep(:,:,cam)*Heuc;
end

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


% ------------------------------------------------------------------------
% 6a. Projective Bundle Adjustment. Use BAProjectiveCalib function
% Coordinates of 3D and 2D points are given in homogeneus coordinates
% ------------------------------------------------------------------------
% auxiliary matrix that indicates that all points are visible in all the cameras

vp = ones(npoints,ncam);
[P_bundle,Q_bundle] = BAProjectiveCalib(q_data,P_rep,Q_2cam_est,vp);

% ------------------------------------------------------------------------
% 6b. Compute the statistics of the reprojection error for the improved projective reconstruction
% ------------------------------------------------------------------------

disp(['Mean re-projection error of initial reconstruction (9 cams)   = ' num2str( ErrorRetroproy(q_data,P_bundle, Q_bundle)/2 )]);
draw_reproj_error(q_data,P_bundle, Q_bundle);

% ------------------------------------------------------------------------
% 7. Proceed as in previous lab to extract the metric reconstruction from
% the optimized projective reconstruction (bundle adjustment)
% ------------------------------------------------------------------------
p_l = zeros(3,6,ncam);
    % we need three pairs of parallel lines
for cam = 1:ncam
    % the first pair, the three faces of the three pairs form a corner in the
    % cube
    p_l(:,1,cam) =cross(q_data(:,1,cam), q_data(:,5,cam));
    p_l(:,2,cam) =cross(q_data(:,2,cam), q_data(:,6,cam));

    % the second pair
    p_l(:,3,cam) =cross(q_data(:,1,cam), q_data(:,2,cam));
    p_l(:,4,cam) =cross(q_data(:,5,cam), q_data(:,6,cam));

    % the third pair
    p_l(:,5,cam) =cross(q_data(:,1,cam), q_data(:,4,cam));
    p_l(:,6,cam) =cross(q_data(:,5,cam), q_data(:,8,cam));
end

v_p = zeros(3,3,ncam); 
for cam =1:ncam
    % we need three vanishing points
    v_p(:,1,cam) = cross(p_l(:,1,cam),p_l(:,2,cam));
    v_p(:,2,cam) = cross(p_l(:,3,cam),p_l(:,4,cam));
    v_p(:,3,cam) = cross(p_l(:,5,cam),p_l(:,6,cam));
    
    %normalize the points by the third coordinate 
    v_p(:,1,cam) =v_p(:,1,cam)/v_p(3,1,cam);
    v_p(:,2,cam) =v_p(:,2,cam)/v_p(3,2,cam);
    v_p(:,3,cam) =v_p(:,3,cam)/v_p(3,3,cam);
end
for k=1:ncam
    figure();
    hold on
    scatter(q(1,:,k),q(2,:,k),30,[1,0,0]);
    draw_proj_cube(q(:,:,k));
    % draw the projection of lines and vanishing points in each image
    scatter(v_p(1,:,k),v_p(2,:,k),'r');
    draw_lines(p_l(:,:,k),-8000,10000,-8000,8000,'b')
    % ...

end

% Obtain the 3D points at the infinity plane by linear triangulation. Use
% the linear_triangulation() function
v_p_3d = linear_triangulation(v_p,P_bundle);

% The following function returns the hom coords. of a plane if X(3,4) contains 3 points of that plane 
plane = NumKernel(v_p_3d.'); 

% Build the P3 transformation H_aff that transforms the projective reconstruction into a affine reconstruction
H_aff = [eye(3,3) zeros(3,1); plane'];

% Apply H_aff to points and projection cameras
Qa = H_aff*Q_bundle;
% Pa(:,:,1) = Q_bundle(:,:,1)*inv(H_aff);
% Pa(:,:,2) = Q_bundle(:,:,2)*inv(H_aff);

draw_3Dcube(Qa);

%Euclidean reconstruction

Q5 = zeros(4,5); % we need 5 points, the first three should be vanishing points, then two opposite vertices for the cube
Q5(:,1:3) = v_p_3d; %4*3 matrix
Q5(:,4)= Q_bundle(:,1); %4*1 coords for the 1st vertice 
Q5(:,5)= Q_bundle(:,7); %4*1 coords for the 7th vertice
Heuc = calc_reference_homography(Q5);

Q_euc = inv(Heuc)*Q_bundle;
for cam = 1:ncam
    P_euc(:,:,cam) = P_bundle(:,:,cam)*Heuc;
end

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

% ------------------------------------------------------------------------
% 8. Extra angle_3points()
% ------------------------------------------------------------------------
