%% Segmentation of Immunofluorescence and Imaging Mass Cytometry nuclei images
%  Author: Nicolas Paolo Canete (ncan4475@uni.sydney.edu.au; nickcee900@gmail.com)
%
%  Modified from the IdentifyPrimAuto module and IdentifySecondary module
%  from CellProfiler (Carpenter et al., 2006).

function segmentedImage = nucleiSegment(nuclei, minDiameter, maxDiameter)
%% Parameters:
%  nuclei: Nuclei stained image to segment.
%
%  minDiameter & maxDiameter: Minimum and maximum diameter of nuclei
%                             expected to be identified in pixels.
%
%  distanceToDilate: Number of pixel dilations to perform to capture cell
%                    body.

%% Usage example
%  nuclei = imread('nuclei.tif');
%  segmentedImage = nucleiSegment(nuclei, 10, 50, 6);

%% Identify Primary Objects
%  Identifies objects given only an image as input. Modified from
%  IdentifyPrimAuto module in CellProfiler.
%
%  Overview of the strategy (adapted from IdentifyPrimAutomatic
%  documentation):
%  Properly identifying primary objects (nuclei) that are well-dispersed,
%  non-confluent, and bright releative to the background is straightforward
%  by applying a simple threshold to the image. This is fast, but usually
%  fails when nuclei are touching. Automatic threshold methods can be used
%  (adaptive Otsu's (Otsu, 1979) is used here). For most biological images,
%  at least some nuclei are touching, so a novel modular three-step
%  strategy from CellProfiler based on previously published algorithms
%  (Malpica et al., 1997; Meyer and Beucher, 1990; Ortiz de Solorzano et
%  al., 1999; Wahlby, 2003; Wahlby et al., 2004) can be used.

%%%%%%%%%%
% Step 1 %
%%%%%%%%%%
%  Objects can be distinguished as an individual nucleus or two or more
%  clumped nuclei. When nuclei are bright in the middle and dimmer towards
%  the edges (the most common case), identifying local maxima in the
%  smoothed intensity image works well.

originalImage = double(nuclei);
threshold = otsuLocal(originalImage);

% Apply a slight smoothing before thresholding to remove 1-pixel objects
% and to smooth the edges of the objects (could be set as user option).
sigma = 1;
filtLength = 2*sigma;
[x,y] = meshgrid(-filtLength:filtLength,-filtLength:filtLength); % Filter kernel grid
f = exp(-(x.^2+y.^2)/(2*sigma^2));
f = f/sum(f(:)); % Gaussian kernel filter
blurredImage = conv2(originalImage,f,'same') ./ conv2(ones(size(originalImage)),f,'same');

objects = blurredImage > threshold;

threshold = mean(threshold(:));

% Fill holes (Note: We do not need to fill holes (might implement as user 
% option).
% objects = imfill(double(objects), 'holes');

% Smooth images to supress maxima (could be set as user option).
sizeOfSmoothingFilter = 2.35*minDiameter/3.5; % FWHM
blurredImage = gaussSmooth(originalImage, sizeOfSmoothingFilter);

% Get local maxima. This will be done on a lower-resolution image for
% speed. The image is resized to a size where the smallest objects are
% about 10 pixels wide. Local maxima within a radius of 5-6 pixels are
% then extracted (could be made into a user-set parameter). This maxima
% supression size should be equal to the minimum acceptable radius if the
% objects are perfectly circular with local maxima in the centre. In
% practice, the min diameter is divided by 1.5 to allow the local maxima to
% be shifted somewhat from the centre of the object.
imageResizeFactor = 10/minDiameter;
maximaSupressionSize = 7; % 10/1.5
maximaMask = getnhood(strel('disk', maximaSupressionSize));
resizedBlurredImage = imresize(blurredImage, imageResizeFactor, 'bilinear');
maximaImage = resizedBlurredImage;
% Extract local maxima
maximaImage(resizedBlurredImage < ordfilt2(resizedBlurredImage, sum(maximaMask(:)), maximaMask)) = 0;
maximaImage = imresize(maximaImage, size(blurredImage), 'bilinear'); % Resize
maximaImage = maximaImage > threshold;
if all(maximaImage(:)) % All > 0 i.e. uniform
    maximaImage = zeros(size(maximaImage));
else
    maximaImage = bwmorph(maximaImage, 'shrink', inf); % Objects shrunk to points
end
% No impose maxima - transforms image so that maxima of objects are found
% where maximaImage is true
overlaid = imimposemin(1-originalImage, maximaImage);

%%%%%%%%%%
% Step 2 %
%%%%%%%%%%
%  The edges of nuclei are identified. For nuclei within the image that do
%  not appear to touch, the edges are easily determined using thresholding.
%  For nuclei that do appear to touch, an intensity approach works. Wheere
%  the dividing lines tend to be dimmer than the remainder of the nucleus
%  (the most comon case), already identified nuclear markers are starting
%  points for a watershed algorithm (Vincent and Soille, 1991) applied to
%  the original image.

% Calculate the watershed transform and cut objects along the boundaries
watershedBoundaries = watershed(overlaid) > 0;
objects = objects.*watershedBoundaries;
objects = bwlabel(objects);

%%%%%%%%%%
% Step 3 %
%%%%%%%%%%
%  Some identified nuclei are discarded or merged together (if required by
%  user). Note that some filters are not required for our analysis and are
%  hence not included here.

% Remove objects with no marker in them (happens occasionally).
pixelIdxList = regionprops(objects, 'PixelIdxList'); %#ok<MRPBW> I don't know why this is here
for i = 1:length(pixelIdxList)
    if sum(maximaImage(pixelIdxList(i).PixelIdxList)) == 0
        objects(pixelIdxList(i).PixelIdxList) = 0;
    end
end

% Get object diameters
diameters = regionprops(objects, 'EquivDiameter');
diameters = [0; cat(1, diameters.EquivDiameter)];

% Remove objects outside of specified diameter range (can be user input).
diameterMap = diameters(objects+1); % Create image with object intensity equal to the diameter.
objects(diameterMap < minDiameter) = 0; % Remove objects that are too small.
prelimPrimLabelMatrix = objects; % For use in Identify Second Objects
objects(diameterMap > maxDiameter) = 0; % Remove objects that are too big.

% Relabel objects
objects = bwlabel(objects > 0);
finalPrimLabelMatrix = objects;

segmentedImage = finalPrimLabelMatrix;
end

%% Otsu's method of thresholding (adapted from CPThreshold subfunction) 
function thresholdVal = otsu(im)
im = im(:);

if max(im) == min(im)
    thresholdVal = im(1);
else
    % Limit dynamic range of the image to 256.
    minVal = max(im)/256;
    im(im < minVal) = minVal;
    im = log(im);
    minVal = min(im);
    maxVal = max(im);
    im = (im - minVal) / (maxVal - minVal);
    thresholdVal = exp(minVal + (maxVal - minVal) * graythresh(im));
end
end

%% Adaptive Otsu's method of thresholding (adapted from CPThreshold subfunction)
function thresholdVal = otsuLocal(im)
linearMaskedImage = im(:);
thresholdVal = otsu(linearMaskedImage);
globalThreshold = thresholdVal;

% Choose the block size that best covers the original image in the sense
% that the number of extra rows and columns is minimal. We want blocks to
% be big enough to contain both background and objects. The minimum block
% size is about 50x50 pixels. The line below divides the image in 10x10
% blocks, and makes sure that the block size is at least 50x50 in pixels.
[m,n] = size(im);
blockSize = max(50,min(round(m/10),round(n/10)));
% Calculates a range of acceptable block sizes as plus-minus 10% of the
% suggested block size.
blockSizeRange = floor(1.1*blockSize):-1:ceil(0.9*blockSize);
[~, minIdx] = min(ceil(m./blockSizeRange).*blockSizeRange - m + ceil(n./blockSizeRange).*blockSizeRange - n);
bestBlockSize = blockSizeRange(minIdx);

% Pads the image so that the blocks fit properly.
rowsToAdd = bestBlockSize*ceil(m/bestBlockSize) - m; % Calculate additional row/columns.
colsToAdd = bestBlockSize*ceil(n/bestBlockSize) - n;
rowsToAddPre = round(rowsToAdd/2); % Calculate how much to add to the start (pre) and end (post).
rowsToAddPost = rowsToAdd - rowsToAddPre;
colsToAddPre = round(colsToAdd/2);
colsToAddPost = colsToAdd - colsToAddPre;
paddedImage = padarray(im, [rowsToAddPre colsToAddPre], 'replicate', 'pre');
paddedImage = padarray(paddedImage, [rowsToAddPost colsToAddPost], 'replicate', 'post');
paddedImageAndCropMask = paddedImage;

% Perform otsu threshold on blocks of size as determined above.
block = [bestBlockSize bestBlockSize];
thresholdVal = blkproc(paddedImageAndCropMask, block, @otsu); %#ok<DBLKPRC> blproc works but not blockproc?

% Resizs the block-produced image to be the size of the padded image.
% Bilinear prevents dipping below zero.
thresholdVal = imresize(thresholdVal, size(paddedImage), 'bilinear');
% The image is cropped to get rid of the padding, to make the result the
% same size as the original image.
thresholdVal = thresholdVal(rowsToAddPre+1:end-rowsToAddPost, colsToAddPre+1:end-colsToAddPost);

% Adjusts any of the threshold values that are significantly lower or
% higher than the global threshold. If there are no objects within a block
% (e.g. if cells are very sparse), an unreasonable threshold will be
% overridden.
thresholdVal(thresholdVal <= 0.7*globalThreshold) = 0.7*globalThreshold;
thresholdVal(thresholdVal >= 1.5*globalThreshold) = 1.5*globalThreshold;
end

%% Gaussian smoothing (adapted from CPSmooth subfunction)
function smoothedImage = gaussSmooth(im, sizeOfSmoothingFilter)
% If the smoothing filter (S) is larger than a predefined effective maximum
% filter size diameter (L), then rescale the original image by L/S and
% rescale S to L.
maxSizeOfSmoothingFilter = 50;
if sizeOfSmoothingFilter >= maxSizeOfSmoothingFilter
    resizingFactor = maxSizeOfSmoothingFilter/sizeOfSmoothingFilter;
    [originalRow,originalCol] = size(im);
    im = imresize(im, resizingFactor);
    sizeOfSmoothingFilter = maxSizeOfSmoothingFilter;
    resized = 1;
else
    resized = 0;
end

% Apply Gaussian filter
sigma = sizeOfSmoothingFilter/2.35; % FWHM to sigma
h = fspecial('gaussian', [round(sizeOfSmoothingFilter) round(sizeOfSmoothingFilter)], sigma); % Generate filter
smoothedImage = imfilter(im, h, 'replicate');

% Restore original size
if resized
    smoothedImage = imresize(smootedImage, [originalRow,originalCol]);
end
end
