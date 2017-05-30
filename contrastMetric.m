function [ contrast ] = contrastMetric( img )
%CONTRASTMETRIC
%LI
% Basic sharpness/contrast metric based on average image gradient
%contrast metric
% See Diatom autofocusing in brightfield microscopy: a comparative study
% Also see A Ringing Metric to Evaluate the Quality of Images Restored using Iterative Deconvolution Algorithms
img = rgb2gray(img);

% Sobel operator
S_x = [1 0 -1; 2 0 -2; 1 0 -1];
S_y = [1 2 1; 0 0 0; -1 -2 -1];
% Convolve with Sobel operator
G_x = imfilter(img,S_x,'conv');
G_y = imfilter(img, S_y, 'conv');
G = sqrt(G_x.^2 + G_y.^2);
meanGrad = mean(G(:));
contrast =(G - meanGrad).^2;
contrast = sum(contrast(:));
end

