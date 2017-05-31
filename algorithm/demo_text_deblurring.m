
% adapted from code by 
% Pan, Jinshan, et al. "Deblurring text images via L0-regularized intensity and gradient prior." Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition. 2014.
%  https://sites.google.com/site/jspanhomepage/l0rigdeblur
% alpha = (0,1) hyperlaplacian
% beta = controls the weighting on the 1st regularization term, bigger more blur, smaller sharper




addpath(genpath('image'));
addpath(genpath('whyte_code'));
addpath(genpath('cho_code'));
opts.prescale = 1; %%downsampling
opts.xk_iter = 5; %% the iterations
opts.gamma_correct = 1.0;
opts.k_thresh = 20;

% UNCOMMENT TO RUN

% test original Pan algorithm run time on computer
%filename = 'image\8_patch_use.png'; opts.kernel_size = 135;  saturation = 0;
%lambda_pixel = 4e-3; lambda_grad = 4e-3; opts.gamma_correct = 2.2;
%lambda_tv = 0.002; lambda_l0 = 2e-4; weight_ring = 1;

% UNCOMMENT TO RUN

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filename = 'image\8_patch_use.png'; opts.kernel_size = 135;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 200; opts.gamma_correct = 2.2;
% lambda_tv = 0.002; lambda_l0 = 2e-4; weight_ring = 1;
% alpha = 0.1; beta = 200;
% 
% filename = 'image\8_patch_use.png'; opts.kernel_size = 135;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 200; opts.gamma_correct = 2.2;
% lambda_tv = 0.002; lambda_l0 = 2e-4; weight_ring = 1;
% alpha = 0.9; beta = 200;
% 
% filename = 'image\8_patch_use.png'; opts.kernel_size = 135;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 200; opts.gamma_correct = 2.2;
% lambda_tv = 0.002; lambda_l0 = 2e-4; weight_ring = 1;
% alpha = 0.1; beta = 300;

% filename = 'image\8_patch_use.png'; opts.kernel_size = 135;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 200; opts.gamma_correct = 2.2;
% lambda_tv = 0.002; lambda_l0 = 2e-4; weight_ring = 1;
% alpha = 0.9; beta = 300;

% filename = 'image\8_patch_use.png'; opts.kernel_size = 135;  saturation = 0;
% lambda_pixel = 4e-3; lambda_grad = 200; opts.gamma_correct = 2.2;
% lambda_tv = 0.002; lambda_l0 = 2e-4; weight_ring = 1;
% alpha = 0.5; beta = 300;

% filename = 'boat2.jpg'; opts.kernel_size = 135;  saturation = 1;
% lambda_pixel = 4e-3; lambda_grad = 200; opts.gamma_correct = 2.2;
% lambda_tv = 0.002; lambda_l0 = 2e-4; weight_ring = 1;
% alpha = 0.1; beta = 200;
% 
% filename = 'boat2.jpg'; opts.kernel_size = 135;  saturation = 1;
% lambda_pixel = 4e-3; lambda_grad = 200; opts.gamma_correct = 2.2;
% lambda_tv = 0.002; lambda_l0 = 2e-4; weight_ring = 1;
% alpha = 1/2; beta = 200;
% 
% filename = 'boat2.jpg'; opts.kernel_size = 135;  saturation = 1;
% lambda_pixel = 4e-3; lambda_grad = 200; opts.gamma_correct = 2.2;
% lambda_tv = 0.002; lambda_l0 = 2e-4; weight_ring = 1;
% alpha = 0.9; beta = 200;


%===================================
y = imread(filename);
% y = y(3:end-2,3:end-2,:);
%y = imfilter(y,fspecial('gaussian',5,2),'same','replicate'); 


if size(y,3)==3
    yg = im2double(rgb2gray(y));
else
    yg = im2double(y);
end
tic;
%%% use different image estimation step
[kernel, interim_latent] = blind_deconv_hyperl(yg, lambda_pixel, lambda_grad, opts, alpha, beta);
%[kernel, interim_latent] = blind_deconv(yg, lambda_pixel, lambda_grad, opts);
toc
y = im2double(y);
%% Final Deblur: 
if ~saturation
    %% 1. TV-L2 denoising method
    Latent = ringing_artifacts_removal(y, kernel, lambda_tv, lambda_l0, weight_ring);
else
    %% 2. Whyte's deconvolution method (For saturated images)
    Latent = whyte_deconv(y, kernel);
end
figure; imshow(Latent)
%%
k = kernel - min(kernel(:));
k = k./max(k(:));
imwrite(k,['results\' filename(7:end-4) '_kernel.png']);
imwrite(Latent,['results\' filename(7:end-4) '_result.png']);
imwrite(interim_latent,['results\' filename(7:end-4) '_interim_result.png']);
%%

