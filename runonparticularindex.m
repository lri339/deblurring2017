function [ score, details ] = runonparticularindex(indices)
%RUNFUNCTIONONIMAGES 
% LI
% not used / for debugging
% calculate on specific indices

% real dataset
baseImageDir = 'C:\Users\lulu\Desktop\Thesis\Develop new metric\cvpr16_deblur_study_all_deblurred_results\real';
deblurredDir = 'C:\Users\lulu\Desktop\Thesis\Develop new metric\cvpr16_deblur_study_real_dataset\real_dataset';

% Get the list of image names
imageNamesFile = fopen('C:\Users\lulu\Desktop\Thesis\Develop new metric\cvpr16_deblur_study-master\cvpr16_deblur_study-master\list\real.txt', 'r');
imageNamesList = textscan(imageNamesFile,'%s %*[^\n]');
fclose(imageNamesFile);
imageNamesList = imageNamesList{1};

% Get the list of method names
methodNamesFile = fopen('C:\Users\lulu\Desktop\Thesis\Develop new metric\cvpr16_deblur_study-master\cvpr16_deblur_study-master\list\method2.txt', 'r');
methodNamesList = textscan(imageNamesFile,'%s %*[^\n]');
fclose(methodNamesFile);
methodNamesList = methodNamesList{1};


% initialise score
score = zeros (1,length(indices));
iter = 1;

for i = 1:length(indices)
    
    %Get the image number
    image_index = floor((indices(i)-1)/14)+1;
    %Get the method number
    method_index = mod(indices(i)-1,14);
    % Concatenate to get image folder name
    currentImageDir = strcat(baseImageDir, '\', imageNamesList(image_index));
    % Check if the folder exists (it should)
    if exist(currentImageDir{1}, 'dir') ~=7
        fprintf('Folder %s does not exist.\n', currentImageDir{1});
    end
    
    %%% RUN FUNCTION ON BLURRY IMAGE
    deblurredFileName = strcat(deblurredDir, '\', imageNamesList(image_index));
    blurred = imread(deblurredFileName{1}, 'jpeg');
    blurred = im2double(blurred);
    
    if method_index==0   
        %try
            [score(i), details(i)] = measure_adapted(blurred,blurred);
        %catch
         %   disp('Error.');
        %end
        iter = iter + 1
    else
        imageFile = strcat(currentImageDir, '\', imageNamesList(image_index), '_', methodNamesList(method_index), '.png');
        if exist(imageFile{1}, 'file')~= 2
            % If image is missing put NaN
            fprintf('File %s is missing.\n', imageFile{1});
            score(iter) = NaN;
        else

            % RUN FUNCTION ON DEBLURRED IMAGE
            deblurred = imread(imageFile{1});
            deblurred = im2double(deblurred);
            %try
                [score(i), details(i)] = measure_adapted(deblurred,blurred);
            %catch
             %   disp('Error.');
            %end
            
        end
        iter = iter + 1
                
    end
end


end



