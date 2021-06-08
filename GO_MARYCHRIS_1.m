% Calibration of cameras
clear all

[checkerboardPoints img_pattern] = get_real_points_checkerboard_vmmc (9, 190, 1);

xyt = zeros ([7, 2, 9]);
H = zeros ([7, 3, 3]);
images = 'D:\Documents\VMMC\FINAL_PROJECT\720\';
figure(1)

hold on

for i = 1:7 
    img = imread (sprintf('%s%d.jpg', images, i));
    xyt(i, :, :) = get_user_points_vmmc (img);
    H(i, :, :) = homography_solve_vmmc (checkerboardPoints', squeeze(xyt(i, :, :)));
    H(i, :, :) = homography_refine_vmmc (checkerboardPoints', squeeze(xyt(i, :, :)), H(i, :, :));
    subplot (7, 2, 2 * i - 1)
    imshow (img)
    tr_image = imtransform (img_pattern, maketform ('projective', squeeze(H(i, :, :))'), 'XData', [1 size(img, 2)], 'YData', [1 size(img, 1)]);  
    subplot (7,2,2*i)
    imshow (tr_image)
end

H2 = {};
for i = 1:7
    H2{i} = squeeze(H(i, :, :));
end
A_small = internal_parameters_solve_vmmc (H2');