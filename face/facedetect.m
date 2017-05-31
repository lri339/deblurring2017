% LI
% run face detection then deblur the face 

%I = imread('Original_255.jpg');
I = imread('motion0034.jpg');
FDetect = vision.CascadeObjectDetector;

 %detect Face
FDetect = vision.CascadeObjectDetector;

%Read the input image
I = imread('motion0034.jpg');

%Returns Bounding Box values based on number of objects
BB = step(FDetect,I);

figure,
imshow(I); hold on

I = imread('motion0034.jpg');
imshow(I);
hold on;
rectangle('Position',BB(2,:),'LineWidth',5,'LineStyle','-','EdgeColor','r');

%%%%%%%%% make region of interest slightly larger than bounding box with face
hold off;
x = BB(2,1)
y = BB(2,2)
w = BB(2,3)
h = BB(2,4)
% w =floor(w/2)
% h = floor(h/2)
%I(y:(y+w),x:(x+h),:) = 1;
imshow(I);
ROI = I(y:(y+w),x:(x+h),:);

[imgDim1, imgDim2, ~] = size(I);

enlargefactorH = 0.2;
enlargefactorV = 0;
y1 = max(1, y-floor(w+enlargefactorH*w));
y2 = min(y+floor(w+w*enlargefactorH), imgDim1);
x1 = max(1, x-floor(h+enlargefactorV*h));
x2 = min(x+floor(h+h*enlargefactorV), imgDim2);

ROI = I(y1:y2,x1:x2,:);
imshow(ROI)
set(gcf,'color','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Run compute_gradient_similarity for the image


%run "Copy_of_computer_gradient_similarity" to find best match
% 002_03_02_051_05.png
% is the bset match
	match_name = '002_03_02_051_05.png'



%%%%%%%%%%%%%%%%%%%%%%%%%%%% Run Deblurring
%%%%%%%%%% runs deblurring
%%%%%%%%% adapted from Pan 2014 code. Need this code. See README.

	y_color = ROI;
	k_size = 30;


      
      blurred = rgb2gray(y_color);
      
      blurred = im2double(blurred);
     
      Matched = imread(['./find_structures_code/Training/' match_name]);
      Matched = imresize(Matched,[size(blurred,1), size(blurred,2)],'bilinear');
      maskname = match_name(1:end-4);
      maskname = [maskname '_mask.png'];
      Mask = imread(['./find_structures_code/Training_mask/' maskname]);
      Mask = imresize(Mask,[size(blurred,1), size(blurred,2)],'bilinear');
      opts.kernel_size  =k_size;
      opts.xk_iter = 50;
      opts.gamma_correct = 1.0;
      
	  %% Deblurring
      [interim_latent, kernel] = blind_deconv(blurred, Matched, Mask, opts);
      
	  %% Pyramid - not used
      %[interim_latent, kernel] = blind_deconv_coarse_to_fine(blurred, Matched, Mask, opts);
	  
	  
	  
	  
	  
	  
	  
      %% Final deblurring
		y_color = double(y_color)/255;
      
      % different tuning
      deblur = [];
      for cc = 1:size(y_color,3)
          deblur(:,:,cc)=deconvSps(y_color(:,:,cc),kernel,0.03,100);
          %deblur(:,:,cc)=L0Restoration(y_color(:,:,cc), kernel, 1000)
      end

		figure(2);
		imshow(deblur);
		set(gcf,'color','w');



		% different tuning

      deblur = [];
      for cc = 1:size(y_color,3)
          deblur(:,:,cc)=deconvSps(y_color(:,:,cc),kernel,0.005,100);
      end

		figure(2);
		imshow(deblur);

		% different tuning


      for cc = 1:size(y_color,3)
          deblur(:,:,cc)=L0Restoration(y_color(:,:,cc), kernel, lambda, 2)
      end

      figure(2);
      imshow(deblur);


















