
%RUNFUNCTIONONIMAGES 
% LI
% get the original LR image blur metric and the related features for the CVR real image dataset 

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
score = zeros (1,1400);
iter = 1;

% Iterate over the datasets for each individual image
for i = 1:length(imageNamesList) 

    % Concatenate to get image folder name
    currentImageDir = strcat(baseImageDir, '\', imageNamesList(i));

    % Check if the folder exists (it should)
    if exist(currentImageDir{1}, 'dir') ~=7
        fprintf('Folder %s does not exist.\n', currentImageDir{1});
    end

    %%% RUN FUNCTION ON BLURRY IMAGE
    deblurredFileName = strcat(deblurredDir, '\', imageNamesList(i));
    blurred = imread(deblurredFileName{1}, 'jpeg');
    blurred = im2double(blurred);
    try
        [score(iter), details(iter)] = measure(blurred,blurred);
    catch
%         score(iter) = NaN;
%         details(iter) = NaN;
    end
    iter = iter + 1

    % iterate over all the method types
    for j = 1:length(methodNamesList)

        imageFile = strcat(currentImageDir, '\', imageNamesList(i), '_', methodNamesList(j), '.png');
        if exist(imageFile{1}, 'file')~= 2
            % If image is missing put NaN
            fprintf('File %s is missing.\n', imageFile{1});
            score(iter) = NaN;
        else

            % RUN FUNCTION ON DEBLURRED IMAGE
            deblurred = imread(imageFile{1});
            deblurred = im2double(deblurred);
            try
                [score(iter), details(iter)] = measure(deblurred,blurred);
            catch
%                 score(iter) = NaN;
%                 details(iter) = NaN;
            end
            %score(iter) = measure(im2double(deblurred), im2double(blurred));
        end
        iter = iter + 1

    end

end





