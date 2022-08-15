function [MatrixLabel,spots_tiff]=  IdentifySpots2D(RNAtiff, iHsize, iImgLimes, ...
    iRescaleThr, iDetectionThr, iDeblendSteps, iObjIntensityThr)
% IDENTIFYPOTS2D: This function is modified from CellProfiler and Battich,
% Stoeger and Pelkmans
% Website: http://www.imls.uzh.ch/research/pelkmans.html
%
% SHORT DESCRIPTION:
% Detects spots as destribed by Battich et al., 2013.
% ***********************************************
% Will Determine Spots in 2D Image stacks after Laplacian Of Gaussian (LoG)
% enhancing of spots. Many of the input arguments are optional. Note that
% while an external script has to be run in order to choose robust values,
% manual selection of the parameters can often yield good estimates, if
% the signal is clear enough.
%
% WHAT DID YOU CALL THE IMAGES YOU WANT TO PROCESS?
% Object detection should be done on this image.
%
% HOW DO YOU WANT TO CALL THE OBJECTS IDENTIFIED PRIOR TO DEBLENDING?
% This is the name of the the spots identified after thresholding the LoG
% image.
%
% HOW DO YOU WANT TO CALL THE OBJECTS IDENTIFIED AFTER DEBLENDING?
% Optional. Deblending can be done after spot detection to separate close
% objects. The algorithm is based upon SourceExtractor. To skip this step,
% insert / as name of the object.
%
% OBJECTSIZE
% This value corresponds to the approximate size of you spots. It should
% be their diameter in pixels. The LoG will use a mask of this size to
% enhance radial signal of that size. Note that in practice the specific value
% does not affect the number of spots, if spots are bright (eg. pixel size 5
% or 6).
%
% INTENSITY QUANTA PER IMAGE
% Prior to spot detection the images are rescaled according to their
% intensity. Since the specific value of minimal and maximal intensities
% are frequently not robust across multiple images, intensity quantile are
% used instead. [0 1] would correspond to using the single dimmest pixel
% for minimal intensity and the single brightest pixel for maximal
% intensity. [0.01 0.90] would mean that the minimum intensity is derived
% from the pixel, which is the 1% brightest pixel of all and that the
% maximum intensity is derived from the pixel, which is the 90% brightest
% pixel .
%
% INTENSITY BORERS FOR INTENSITY RESCALING OF IMAGES
% Most extreme values that the image intensity minimum and image intensity
% maximum (as defined by the quanta) are allowed to have
% [LowestPossibleGreyscaleValueForImageMinimum
% HighestPossibleGreyscaleValueForImageMinimum
% LowestPossibleGreyscaleValueForImageMaximum
% HighestPossibleGreyscaleValueForImageMaximum]
% To ignore individual values, place a NaN.
% Note that these parameters very strongly depend upon the variability of
% your illumination source. When using a robust confocal microscope you can
% set the lowest and highest possible values to values,  which are very
% close (or even identical). If your light source is variable during the
% acquisition (which can be the case with Halogen lamps) you might choose
% less strict borders to detect spots of varying intensites.
%
% THRESHOLD OF SPOT DETECTION
% This is the threshold value for spot detection. The higher it is the more
% stringent your spot detection is. Use external script to determine a
% threshold where the spot number is robust against small variations in the
% threshold.
%
% HOW MANY STEPS OF DEBLENDING DO YOU WANT TO DO?
% The amount of deblending steps, which are done. The higher it is the less
% likely it is that two adjacent spots are not separated. The default of 30
% works very well (and we did not see improvement on our images with higher
% values). Note that the number of deblending steps is the main determinant
% of computational time for this module.
%
% WHAT IS THE MINIMAL INTENSITY OF A PIXEL WITHIN A SPOT?
% Minimal greyscale value of a pixel, which a pixel has to have in order to
% be recognized to be within a spot. Opitonal argument to make spot
% detection even more robust against very dim spots. In practice, we have
% never observed that this parameter would have any influence on the spot
% detection. However, you might include it as an additional safety measure.
%
% WHICH IMAGE DO YOU WANT TO USE AS A REFERENCE FOR SPOT BIAS CORRECTION?
% Here you can name a correction matrix which counteracts bias of the spot
% correction across the field of view. Note that such a correction matrix
% has to be loaded previously by a separate module, such as
% LOADSINGLEMATRIX
%
% The design of this module largely follows a IdentifyPrimLoG2 by
% Baris Sumengen.
%
% $Revision: 1889 $



% Intensity Quanta Of Image
[isSafe iImgLimes]= inputVectorsForEvalCP3D(iImgLimes,true);

% Rescale Thresholds
[isSafe iRescaleThr]= inputVectorsForEvalCP3D(iRescaleThr,true);

%As long as user isn't satisfied with parameters keep looping
% Initiate Settings for deblending
Options.ObSize = iHsize;
Options.limQuant = eval(iImgLimes);
Options.RescaleThr = eval(iRescaleThr);
Options.ObjIntensityThr = [];
Options.closeHoles = false;
Options.ObjSizeThr = [];
Options.ObjThr = iDetectionThr;
Options.StepNumber = iDeblendSteps;
Options.numRatio = 0.20;
Options.doLog = 0;


%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGE ANALYSIS %%%
%%%%%%%%%%%%%%%%%%%%%%

currTiff = RNAtiff;
Image = double(currTiff);%.*65535; %convert to scale used for spotdetection
op = fspecialCP3D('2D LoG',iHsize); % force 2D filter



% Object intensity Threshold
% if iObjIntensityThr == '/'
%     iObjIntensityThr = [];
% elseif strcmp(iObjIntensityThr,'auto')
%     iObjIntensityThr = multithresh(Image).*0.7;% be sensitive enough!
% else
%     try
%         iObjIntensityThr = str2double(iObjIntensityThr);
%     catch errObjIntensityThr
%         error(['Image processing was canceled in the ', ModuleName, ' module because the Stepsize for deblending could not be converted to a number.'])
%     end
% end

iObjIntensityThr = multithresh(Image).*0.7;

DetectionBias = [];
% Detect objects, note that input vectors are eval'ed
[ObjCount{1} SegmentationCC{1} FiltImage] = ObjByFilter(Image,op,iDetectionThr,eval(iImgLimes),eval(iRescaleThr),iObjIntensityThr,true,[],DetectionBias);
% Convert to CP1 standard: labelmatrix
MatrixLabel{1} = double(labelmatrix(SegmentationCC{1}));

% Security check, if conversion is correct
if max(MatrixLabel{1}(:)) ~= ObjCount{1}
    error(['Image processing was canceled in the ', ModuleName, ' module because conversion of format of segmentation was wrong. Contact Thomas.'])
end

% Deblend objects
if iDeblendSteps > 0        % Only do deblending, if number of iterations was defined
    MatrixLabel{2} = SourceExtractorDeblend(Image,SegmentationCC{1},FiltImage,Options);
    ObjCount{2} = max(MatrixLabel{2}(:));
end

try
    spots_tiff{1} = MatrixLabel{2} ~= 0;
catch
    spots_tiff{1} = MatrixLabel{1} ~= 0;
end

