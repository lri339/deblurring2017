function [eout,thresh,gv_45,gh_135] = cannyModified(varargin)
%%
%% This is MATLAB's default canny edge finder modified to use bilateral filter to preprocess 
% derivatives instead of gaussian


[a,method,thresh,sigma,thinning,H,kx,ky] = parse_inputs(varargin{:});

% Check that the user specified a valid number of output arguments
if ~any(strcmp(method,{'sobel','roberts','prewitt'})) && (nargout>2)
    error(message('images:edge:tooManyOutputs'))
end

% Transform to a double precision intensity image if necessary
if ~isa(a,'double') && ~isa(a,'single')
    a = im2single(a);
end

[m,n] = size(a);

if strcmp(method,'canny')
    % Magic numbers
    PercentOfPixelsNotEdges = .7; % Used for selecting thresholds
    ThresholdRatio = .4;          % Low thresh is this fraction of the high.
    
    % Calculate gradients using a derivative of Gaussian filter
    [dx, dy] = smoothGradient(a, sigma);
    
    % Calculate Magnitude of Gradient
    magGrad = hypot(dx, dy);
    
    % Normalize for threshold selection
    magmax = max(magGrad(:));
    if magmax > 0
        magGrad = magGrad / magmax;
    end
    
    % Determine Hysteresis Thresholds
    [lowThresh, highThresh] = selectThresholds(thresh, magGrad, PercentOfPixelsNotEdges, ThresholdRatio, mfilename);
    
    % Perform Non-Maximum Suppression Thining and Hysteresis Thresholding of Edge
    % Strength
    e = thinAndThreshold(dx, dy, magGrad, lowThresh, highThresh);
    thresh = [lowThresh highThresh];
    
elseif strcmp(method,'canny_old')
    % Magic numbers
    GaussianDieOff = .0001;
    PercentOfPixelsNotEdges = .7; % Used for selecting thresholds
    ThresholdRatio = .4;          % Low thresh is this fraction of the high.
    
    % Design the filters - a gaussian and its derivative
    
    pw = 1:30; % possible widths
    ssq = sigma^2;
    width = find(exp(-(pw.*pw)/(2*ssq))>GaussianDieOff,1,'last');
    if isempty(width)
        width = 1;  % the user entered a really small sigma
    end
    
    t = (-width:width);
    gau = exp(-(t.*t)/(2*ssq))/(2*pi*ssq);     % the gaussian 1D filter
    
    % Find the directional derivative of 2D Gaussian (along X-axis)
    % Since the result is symmetric along X, we can get the derivative along
    % Y-axis simply by transposing the result for X direction.
    [x,y]=meshgrid(-width:width,-width:width);
    dgau2D=-x.*exp(-(x.*x+y.*y)/(2*ssq))/(pi*ssq);
    
    % Convolve the filters with the image in each direction
    % The canny edge detector first requires convolution with
    % 2D Gaussian, and then with the derivative of a Gaussian.
    % Since Gaussian filter is separable, for smoothing, we can use
    % two 1D convolutions in order to achieve the effect of convolving
    % with 2D Gaussian.  We convolve along rows and then columns.
    
    %smooth the image out
    aSmooth=imfilter(a,gau,'conv','replicate');   % run the filter across rows
    aSmooth=imfilter(aSmooth,gau','conv','replicate'); % and then across columns
    
    %apply directional derivatives
    ax = imfilter(aSmooth, dgau2D, 'conv','replicate');
    ay = imfilter(aSmooth, dgau2D', 'conv','replicate');
    
    mag = sqrt((ax.*ax) + (ay.*ay));
    magmax = max(mag(:));
    if magmax>0
        mag = mag / magmax;   % normalize
    end
    
    % Select the thresholds
    if isempty(thresh)
        counts=imhist(mag, 64);
        highThresh = find(cumsum(counts) > PercentOfPixelsNotEdges*m*n,...
            1,'first') / 64;
        lowThresh = ThresholdRatio*highThresh;
        thresh = [lowThresh highThresh];
    elseif length(thresh)==1
        highThresh = thresh;
        if thresh>=1
            error(message('images:edge:thresholdMustBeLessThanOne'))
        end
        lowThresh = ThresholdRatio*thresh;
        thresh = [lowThresh highThresh];
    elseif length(thresh)==2
        lowThresh = thresh(1);
        highThresh = thresh(2);
        if (lowThresh >= highThresh) || (highThresh >= 1)
            error(message('images:edge:thresholdOutOfRange'))
        end
    end
    
    % The next step is to do the non-maximum suppression. We will accrue
    % indices which specify ON pixels in strong edgemap The array e will become
    % the weak edge map.
    e = cannyFindLocalMaxima(ax,ay,mag,lowThresh);
    
    if ~isempty(e)
        [rstrong,cstrong] = find(mag>highThresh & e);
        
        if ~isempty(rstrong) % result is all zeros if idxStrong is empty
            e = bwselect(e, cstrong, rstrong, 8);
            e = bwmorph(e, 'thin', 1); % Thin double (or triple) pixel wide contours
        end
    end
    
elseif any(strcmp(method, {'log','zerocross'}))
    % The output edge map:
    e = false(m,n);
    
    rr = 2:m-1; cc=2:n-1;
    
    % We don't use image blocks here
    if isempty(H),
        fsize = ceil(sigma*3) * 2 + 1;  % choose an odd fsize > 6*sigma;
        op = fspecial('log',fsize,sigma);
    else
        op = H;
    end
    
    op = op - sum(op(:))/numel(op); % make the op to sum to zero
    b = imfilter(a,op,'replicate');
    
    if isempty(thresh)
        thresh = .75*mean2(abs(b));
    end
    
    % Look for the zero crossings:  +-, -+ and their transposes
    % We arbitrarily choose the edge to be the negative point
    [rx,cx] = find( b(rr,cc) < 0 & b(rr,cc+1) > 0 ...
        & abs( b(rr,cc)-b(rr,cc+1) ) > thresh );   % [- +]
    e((rx+1) + cx*m) = 1;
    [rx,cx] = find( b(rr,cc-1) > 0 & b(rr,cc) < 0 ...
        & abs( b(rr,cc-1)-b(rr,cc) ) > thresh );   % [+ -]
    e((rx+1) + cx*m) = 1;
    [rx,cx] = find( b(rr,cc) < 0 & b(rr+1,cc) > 0 ...
        & abs( b(rr,cc)-b(rr+1,cc) ) > thresh);   % [- +]'
    e((rx+1) + cx*m) = 1;
    [rx,cx] = find( b(rr-1,cc) > 0 & b(rr,cc) < 0 ...
        & abs( b(rr-1,cc)-b(rr,cc) ) > thresh);   % [+ -]'
    e((rx+1) + cx*m) = 1;
    
    % Most likely this covers all of the cases.   Just check to see if there
    % are any points where the LoG was precisely zero:
    [rz,cz] = find( b(rr,cc)==0 );
    if ~isempty(rz)
        % Look for the zero crossings: +0-, -0+ and their transposes
        % The edge lies on the Zero point
        zero = (rz+1) + cz*m;   % Linear index for zero points
        zz = (b(zero-1) < 0 & b(zero+1) > 0 ...
            & abs( b(zero-1)-b(zero+1) ) > 2*thresh);     % [- 0 +]'
        e(zero(zz)) = 1;
        zz = (b(zero-1) > 0 & b(zero+1) < 0 ...
            & abs( b(zero-1)-b(zero+1) ) > 2*thresh);     % [+ 0 -]'
        e(zero(zz)) = 1;
        zz = (b(zero-m) < 0 & b(zero+m) > 0 ...
            & abs( b(zero-m)-b(zero+m) ) > 2*thresh);     % [- 0 +]
        e(zero(zz)) = 1;
        zz = (b(zero-m) > 0 & b(zero+m) < 0 ...
            & abs( b(zero-m)-b(zero+m) ) > 2*thresh);     % [+ 0 -]
        e(zero(zz)) = 1;
    end
 end   
if nargout==0,
    imshow(e);
else
    eout = e;
end







function [I,Method,Thresh,Sigma,Thinning,H,kx,ky] = parse_inputs(varargin)
% OUTPUTS:
%   I      Image Data
%   Method Edge detection method
%   Thresh Threshold value
%   Sigma  standard deviation of Gaussian
%   H      Filter for Zero-crossing detection
%   kx,ky  From Directionality vector

narginchk(1,5)

I = varargin{1};

validateattributes(I,{'numeric','logical'},{'real','nonsparse','2d'},mfilename,'I',1);

% Defaults
Method    = 'sobel';
Direction = 'both';
Thinning  = true;

methods    = {'canny','canny_old','prewitt','sobel','marr-hildreth','log','roberts','zerocross'};
directions = {'both','horizontal','vertical'};
options    = {'thinning','nothinning'};

% Now parse the nargin-1 remaining input arguments

% First get the strings - we do this because the interpretation of the
% rest of the arguments will depend on the method.
nonstr = [];   % ordered indices of non-string arguments
for i = 2:nargin
    if ischar(varargin{i})
        str = lower(varargin{i});
        j = find(strcmp(str,methods));
        k = find(strcmp(str,directions));
        l = find(strcmp(str,options));
        if ~isempty(j)
            Method = methods{j(1)};
            if strcmp(Method,'marr-hildreth')
                error(message('images:removed:syntax','EDGE(I,''marr-hildreth'',...)','EDGE(I,''log'',...)')) 
            end
        elseif ~isempty(k)
            Direction = directions{k(1)};
        elseif ~isempty(l)
            if strcmp(options{l(1)},'thinning')
                Thinning = true;
            else
                Thinning = false;
            end
        else
            error(message('images:edge:invalidInputString', varargin{ i }))
        end
    else
        nonstr = [nonstr i]; %#ok<AGROW>
    end
end

% Now get the rest of the arguments
[Thresh,Sigma,H,kx,ky] = images.internal.parseNonStringInputsEdge(varargin,Method,Direction,nonstr);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Local Function : smoothGradient
%
function [GX, GY] = smoothGradient(I, sigma)

% Create an even-length 1-D separable Derivative of Gaussian filter

% Determine filter length
filterLength = 8*ceil(sigma);
n = (filterLength - 1)/2;
x = -n:n;

% Create 1-D Gaussian Kernel
c = 1/(sqrt(2*pi)*sigma);
gaussKernel = c * exp(-(x.^2)/(2*sigma^2));

% Normalize to ensure kernel sums to one
gaussKernel = gaussKernel/sum(gaussKernel);

% Create 1-D Derivative of Gaussian Kernel
derivGaussKernel = gradient(gaussKernel);

% Normalize to ensure kernel sums to zero
negVals = derivGaussKernel < 0;
posVals = derivGaussKernel > 0;
derivGaussKernel(posVals) = derivGaussKernel(posVals)/sum(derivGaussKernel(posVals));
derivGaussKernel(negVals) = derivGaussKernel(negVals)/abs(sum(derivGaussKernel(negVals)));

% Compute smoothed numerical gradient of image I along x (horizontal)
% direction. GX corresponds to dG/dx, where G is the Gaussian Smoothed
% version of image I.
%filteredIm = bfilter2(I, 6, [3,100/255]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% MODIFICATION IS HERE
%%%%%%
[filteredIm,~] = fastBilateral(im2uint8(I), 100, 3,0.01 , 'Gauss');
filteredIm = double(filteredIm)/255;


GX = imfilter(filteredIm, derivGaussKernel, 'conv', 'replicate');
%GX = imfilter(I, gaussKernel', 'conv', 'replicate'); 
%GX = imfilter(GX, derivGaussKernel, 'conv', 'replicate');

% Compute smoothed numerical gradient of image I along y (vertical)
% direction. GY corresponds to dG/dy, where G is the Gaussian Smoothed
% version of image I.
GY  = imfilter(filteredIm, derivGaussKernel', 'conv', 'replicate');
%GY = imfilter(I, gaussKernel, 'conv', 'replicate');
%GY  = imfilter(GY, derivGaussKernel', 'conv', 'replicate');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Local Function : selectThresholds
%
function [lowThresh, highThresh] = selectThresholds(thresh, magGrad, PercentOfPixelsNotEdges, ThresholdRatio, ~)

[m,n] = size(magGrad);

% Select the thresholds
if isempty(thresh)
    counts=imhist(magGrad, 64);
    highThresh = find(cumsum(counts) > PercentOfPixelsNotEdges*m*n,...
        1,'first') / 64;
    lowThresh = ThresholdRatio*highThresh;
elseif length(thresh)==1
    highThresh = thresh;
    if thresh>=1
        error(message('images:edge:thresholdMustBeLessThanOne'))
    end
    lowThresh = ThresholdRatio*thresh;
elseif length(thresh)==2
    lowThresh = thresh(1);
    highThresh = thresh(2);
    if (lowThresh >= highThresh) || (highThresh >= 1)
        error(message('images:edge:thresholdOutOfRange'))
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Local Function : thinAndThreshold
%
function H = thinAndThreshold(dx, dy, magGrad, lowThresh, highThresh)
% Perform Non-Maximum Suppression Thining and Hysteresis Thresholding of
% Edge Strength
    
% We will accrue indices which specify ON pixels in strong edgemap
% The array e will become the weak edge map.

E = cannyFindLocalMaxima(dx,dy,magGrad,lowThresh);

if ~isempty(E)
    [rstrong,cstrong] = find(magGrad>highThresh & E);
    
    if ~isempty(rstrong) % result is all zeros if idxStrong is empty
        H = bwselect(E, cstrong, rstrong, 8);
    else
        H = false(size(E));
    end
else
    H = false(size(E));
end