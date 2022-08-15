%% Spot Counting Parameters
%ObjectSize:
iHsize = 6;

%Intensity Quanta Per Image:
iImgLimes = '[0.01 0.99]';

%Intensity borders for intensity rescaling of images
%[MinOfMinintens MaxOfMinintens MinOfMaxintens MaxOfMaxintens]
iRescaleThr = '[NaN 120 500 NaN]';

%Threshold of Spot Detection
iDetectionThr = 0.01;

%How many Steps of Deblending do you want to do?
iDeblendSteps = 2;

%What is the minimal intensity of a pixel within a spot?
iObjIntensityThr = NaN;

%% Run
hiv = imread('HIVRNA.tif');

[objects, spots] = IdentifySpots2D(hiv, iHsize, iImgLimes, ...
    iRescaleThr, iDetectionThr, iDeblendSteps, iObjIntensityThr);

%objects = objects{1};
spots = spots{1};

%imwrite(uint16(objects), 'objects.tif');
imwrite(uint8(spots)*255, 'spots.tif');

%Gets row and column where non 0 values are in mask
[spotY, spotX] = find(spots);
csvwrite('HIV.csv',[spotX spotY]);