function Q = darkChannelAverage( img,img_blurry, N, threshold )

% LI
% calculates a score based on dark channel pixel intensity
% higher = bad, lower = better
% img the image: must be colour
% blurryimg
% N is patch isze in pixels
% uses RGB 0 to 255 (8 bit) image representation
% 0 = black
% threshold is percentage over the average minimum



dims = size(img_blurry);
H = dims(1);
W = dims(2);
w = floor(W/N);
h = floor(H/N);
Q = 0;

% find average minimum of blurry image and use to determine threshold
avg_minimum = 0;
threshold = -0.1;

% for m = 1:h
%     for n = 1:w
%         
%         %select patch of image
%         AOI = img_blurry(N*(m-1)+1:N*m, N*(n-1)+1:N*n,:);
%         
%         %find dark cahnnel pixel value
%         darkchannel = min(AOI, [], 3); 
%         
%         avg_minimum = avg_minimum + mean(darkchannel(:));
%         
%     end
% end
% avg_minimum = avg_minimum/(w*h);
% threshold = (1+threshold)*avg_minimum;





for m = 1:h
    for n = 1:w
        
        
        AOI = img_blurry(N*(m-1)+1:N*m, N*(n-1)+1:N*n,:);
        darkchannel = min(AOI, [], 3); 
        threshold = 0.95*median(darkchannel(:));
        
        
        %select patch of image
        AOI = img(N*(m-1)+1:N*m, N*(n-1)+1:N*n,:);
        
        %find dark cahnnel pixel value
        darkchannel = min(AOI, [], 3);
        test = darkchannel > threshold;
        
        % average or sum
        % will try using L0 (i.e. the count of non-zero darkchannel pixels)
        % rather than L1 (i.e. the mean value of darkchannel pixels)
        
        
        %Q = Q + median(darkchannel(:));
        Q = Q + sum(test(:));
        
    end
end







Q = Q/(w*h)/double(N)^2;
end