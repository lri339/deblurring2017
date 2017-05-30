function [ ringingMetric ] = ringingFeatureSimple( img, blurryImg)
%RINGINGFEATURESIMPLE 
% calcualtes simpel ringing metric
%see A Ringing Metric to Evaluate the Quality of Images Restored using Iterative Deconvolution Algorithms

%% Create Mask for Area ringing might be in by finding
%% edges in blurry image
[~, ~, sizeImg] = size(img);
if sizeImg == 3
    img = rgb2gray(img);
    blurryImg = rgb2gray(blurryImg);
end


E_ref = edge(blurryImg, 'canny');
structElement = strel('square',8);
E_mask = imdilate(E_ref, structElement);

%Find edges in deblurred image
E_restored = edge(img, 'canny');

%difference 
tmp1 = E_restored(find(E_mask));
tmp1 = sum(tmp1(:));
tmp2 = sum(E_ref(:));
ringingMetric = (tmp1 - tmp2)/ tmp2 ;


end

