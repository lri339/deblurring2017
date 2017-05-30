function [ ringingMetric ] = ringingFeature( blurred, restored )
%RINGINGFEATURE 
% LI
% blurred = original blurred image
% restored = deblurred image 
% this function calculates a perceptually relevant ringing metric 


img = blurred;


%%%%%%%%%%%% First part: EDGE EXTRACTION %%%%%%%%%%%%%%%%%%
% For speed select grayscale (luminance channel) only.
[~, ~, sizeImg] = size(img);
if sizeImg == 3
    img = rgb2gray(img);
end

[~, ~, sizeImg] = size(restored);
if sizeImg == 3
    restored = rgb2gray(restored);
end

restored = im2double(restored);
img = im2double(img);

% bilateral filter - not used, included within cannyModified()
%img = im2single(img);
% sigmaD = 3;
% sigmaR = 100;
% %spatial standard deviation is 3, therefore choose half window of 3*2
% %this is 95%+ of data.
% halfWindowWidth = 2*sigmaD;
% edgesImg = bfilter2(img, halfWindowWidth, [sigmaD, sigmaR]);

%%%% EDGE FINDING

%canny_thr_hi = 0.85;
canny_thr_hi = 0.75;
% canny_thr_low = 0.4*canny_thr_hi;

% canny is matlab's canny modified to use bilateral smoothing on image
% may not use due to slower speed.
edges = cannyModified(img,'Canny',canny_thr_hi);


%%%%% NOISE REMOVAL %%%%%
%or use alternative code
edges = skeletonize(edges);
% Peter Kovesi function edgelink()
[edgelist edgeim, ~] = edgelink(edges);
% remove elements below a certain number of pixels (use edgelist))
minPixelCount = 20;
i = length(edgelist);
while i>0
    [pixelCount, ~] = size(edgelist{i});
    if pixelCount < minPixelCount
        zeroIndices = sub2ind(size(img), edgelist{i}(:,1), edgelist{i}(:,2));
        edgeim(zeroIndices) = 0;
        edgelist(i)=[];
    end
    i=i-1;
end

%%%%%%%%%%%% Second part: RINGING REGION DETECTION %%%%%%%%%%%%%%%%%%

%make regions
%note: since the matrix is very sparse, this can be optimised by making
%a custom sparse image dilation function
% however, using imdilate for now

%use default 1 pixel edge for edge Region
structElement1 = strel('square',2);
edgeRegion = imdilate(edgeim, structElement1);
%single-sided
%'single-sided support dimension of 4 pixels' = 8 pixel dilated area (?)
% they used 8
structElement2 = strel('square',12);
structElement3 = strel('square',14);
detectionRegion = imdilate(edgeim, structElement2);
% Use same feature and detection region size.
featureRegion = detectionRegion;
%structElement3 = strel('square',17);
%featureRegion = imdilate(edgeim, structElement1);


% Texture Filtering
%Sobel operator on fexreg to find image regions of low texture where
%ringing will be more visible

% why use this? why not just use already computed matrix edges???
%LA = imgradient(img.*(featureRegion>0), 'Sobel');
LA = imgradient(img, 'Sobel');
% range [0.6 0.95]
%textureThreshold = 0.9;
textureThreshold = 0.95;

textureRegion = detectionRegion.*(LA>=textureThreshold).*(~edgeRegion);

% not used - manual implementation
% for k = 1:length(edgelist)
%     i = edgelist{k}(1,:);
%     j = edgelist{k}(2,:);
%     
%     
% end
% for i j in fexreg
% 	la(i,j) = abs(im(i-1,j-1) + 2*im(i-1,j) + im(i-1, j+1)- (im(i+1,j-1)+2*im(i+1,j) + im(i+1,j+1)) )
% 	la(i,j) = la(i,j) + abs(im(i-1, j+1)+2*im(i,j+1)+im(i+1,j+1) - (im(i-1,j-1) + 2*im(i,j-1)+im(i+1,j-1)))
% 	if la(i,j) < textureThreshold
% 		lbm(i,j) = 0 
% 	else
% 		lbm(i,j) = 1;
% 	end
% end

try
    %% clean up: dilate and remove noise.
    textureRegion = imdilate(textureRegion, structElement2);
    [edgelistD, textureRegion, ~] = edgelink(textureRegion);
    % remove elements below a certain number of pixels (use edgelist))
    minPixelCount = 20;
    i = length(edgelistD);
    while i > 0
        [pixelCount, ~] = size(edgelistD{i});
        if pixelCount < minPixelCount
            zeroIndices = sub2ind(size(detectionRegion), edgelistD{i}(:,1), edgelistD{i}(:,2));
            textureRegion(zeroIndices) = 0;
            edgelistD(i)=[];
        end
        i = i-1;
    end

    textureRegion = imdilate(textureRegion, structElement3);
    detectionRegion = (detectionRegion>0) & ~textureRegion;
catch
    % if e.g. textureRegion is empty edgelink() will error
end
% %Luminance filtering
% % ringing not as noticeable in veyr dark or bright regions
% % luminance threshold, between 0 and 0.8
% %ThresholdLum =0.75;
ThresholdLum =0.1;
lowLumiTrue = 0;

for i = 1:length(edgelist)
    
    [sizeList, ~] = size(edgelist{i});
    for j = 1:sizeList
        % local mean luminance
        xCoord = edgelist{i}(j,1);
        yCoord = edgelist{i}(j,2);
        LML = img(xCoord-1:xCoord+1, yCoord-1:yCoord+1);
        LML = sum(LML(:))/9;
        if LML <= ThresholdLum
            edgeim(xCoord, yCoord) = 1;
            lowLumiTrue = 1;
        else
            edgeim(xCoord, yCoord) = 0;
        end
    end
    
end


if lowLumiTrue == 1
    
    %% clean up: dilate and remove noise.
    luminanceRegion = imdilate(edgeim, structElement2);
    [edgelistD, luminanceRegion , ~] = edgelink(luminanceRegion );
    % remove elements below a certain number of pixels (use edgelist))
    minPixelCount = 20;
    i = length(edgelistD);
    while i > 0
        [pixelCount, ~] = size(edgelistD{i});
        if pixelCount < minPixelCount
            zeroIndices = sub2ind(size(detectionRegion), edgelistD{i}(:,1), edgelistD{i}(:,2));
            luminanceRegion (zeroIndices) = 0;
            edgelistD(i)=[];
        end
        i = i-1;
    end

    luminanceRegion = imdilate(luminanceRegion, structElement3);
    detectionRegion = (detectionRegion>0) & ~luminanceRegion;
end



%
% %% clean up: dilate and remove noise.
% detectionRegion = imdilate(detectionRegion, structElement1);
% [edgelistD, ~, ~] = edgelink(detectionRegion);
% % remove elements below a certain number of pixels (use edgelist))
% minPixelCount = 20;
% i = length(edgelistD);
% while i > 0
%     [pixelCount, ~] = size(edgelistD{i});
%     if pixelCount < minPixelCount
%         zeroIndices = sub2ind(size(detectionRegion), edgelistD{i}(:,1), edgelistD{i}(:,2));
%         detectionRegion(zeroIndices) = 0;
%         edgelistD(i)=[];
%     end
%     i = i-1;
% end

%no 'unimpaired regions' removal
% N/A for blurry image

% To determine metric
% take edges selected in detection region
% base this off blurry image for consistnecy (although not ideal)
% and compare the gradient difference:



%Find edges in deblurred image
%E_restored = edge(restored, 'canny');

% this is quite slow due to use of extremely slow bilateral filter
E_restored = cannyModified(restored,'Canny',0.05);
E_restored = skeletonize(E_restored);
%imshow(E_restored)
E_blurry = cannyModified(img,'Canny',0.05);
E_blurry = skeletonize(E_blurry);

%difference in gradients..
% could possibly take max(diff, 0) as in LR paper.
% tmp1 = E_restored(find(detectionRegion));
% tmp1 = sum(tmp1(:));
% tmp2 = edges(find(detectionRegion));
% tmp2 = sum(tmp2(:));
% ringingMetric = (tmp1 - tmp2)/ tmp2 ;

%make sure original edge region left out
detectionRegion = detectionRegion & (~edgeRegion);

ringingIndices = find(detectionRegion);
%tmp1 = edges(ringingIndices);
tmp1 = E_blurry(ringingIndices);
tmp2 = max(0,E_restored(ringingIndices) - tmp1 );
ringingMetric = sum(tmp2(:))/sum(tmp1(:))

% if no strong edges detected, metric is not applicable
if ringingMetric == Inf
    ringingMetric = NaN;
end

% if no detection region detected, set to 0
if isempty(ringingIndices)
    ringingMetric = 0;
end


end

