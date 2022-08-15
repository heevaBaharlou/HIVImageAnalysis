workDir = pwd;
files = dir;
directoryNames = {files([files.isdir]).name};
directoryNames = directoryNames(~ismember(directoryNames,{'.','..'}));

[num,txt,raw] = xlsread('thresholds.csv');

%Manually change values in spreadsheet
% num(32,2) = 2500;
% num(35,2) = 1600;

%Get dimension of file name column
[r,~] = size(txt);
fileNames = cell(1,r-1);

%Store file names in fileNames
for i = 2:r
    fileNames{i-1} = txt{i,1};
end

% toRun = [8];

for i = 1:length(directoryNames)
    disp(i)
    disp(directoryNames{i})
%     try
        %% Filter small nuclear objects
        nuclearMask = logical(imread([workDir filesep directoryNames{i} filesep 'nuclearMask.tif']));
        nuclearMask = bwpropfilt(nuclearMask,'EquivDiameter',[10 50]);
        nuclearMaskLab = bwlabel(nuclearMask);
        
        %% Classify Nuclei
        cd3 = imread([workDir filesep directoryNames{i} filesep 'CD3_R.tif']);
        cd11c = imread([workDir filesep directoryNames{i} filesep 'CD11c.tif']);
        fxiiia = imread([workDir filesep directoryNames{i} filesep 'FXIIIA.tif']);
        
        cd3Blur = imgaussfilt(cd3,1.5);
        cd11cBlur = imgaussfilt(cd11c,1.5);
        fxiiiaBlur = imgaussfilt(fxiiia,1.5);
        
        cd3Thresh = 0;
        cd11cThresh = 0;
        fxiiiaThresh = 0;
        
        for j = 1:length(fileNames)
            if strcmp(fileNames{j},directoryNames{i})
                cd3Thresh = num(j,4);
                cd11cThresh = num(j,3);
                fxiiiaThresh = num(j,2);
            end
        end
        
        cd3Mask = cd3Blur>cd3Thresh;
        cd11cMask = cd11cBlur>cd11cThresh;
        fxiiiaMask = fxiiiaBlur>fxiiiaThresh;
        
        cd11cAF = logical(imread([workDir filesep directoryNames{i} filesep 'CD11c_AFidE.tif']));
        fxiiiaAF = logical(imread([workDir filesep directoryNames{i} filesep 'FXIIIA_AFidE.tif']));
        
        cd11cMask = cd11cMask-cd11cAF;
        fxiiiaMask = fxiiiaMask-fxiiiaAF;
        
        imwrite(uint8(cd3Mask*255), [workDir filesep directoryNames{i} filesep 'cd3Mask.tif']);
        imwrite(uint8(cd11cMask*255), [workDir filesep directoryNames{i} filesep 'cd11cMask.tif']);
        imwrite(uint8(fxiiiaMask*255), [workDir filesep directoryNames{i} filesep 'fxiiiaMask.tif']);
        
%         cd3Val = cell2mat(struct2cell(regionprops(nuclearMask,cd3,'MeanIntensity')));
%         cd11cVal = cell2mat(struct2cell(regionprops(nuclearMask,cd11c,'MeanIntensity')));
%         fxiiiaVal = cell2mat(struct2cell(regionprops(nuclearMask,fxiiia,'MeanIntensity')));
        
        cd3Percentage = cell2mat(struct2cell(regionprops(nuclearMask,cd3Mask,'MeanIntensity'))).*100;
        cd11cPercentage = cell2mat(struct2cell(regionprops(nuclearMask,cd11cMask,'MeanIntensity'))).*100;
        fxiiiaPercentage = cell2mat(struct2cell(regionprops(nuclearMask,fxiiiaMask,'MeanIntensity'))).*100;
        
        allPercentage = [cd3Percentage; cd11cPercentage; fxiiiaPercentage-20];
        %Retrieve the max value in each column (allMax) and the row number
        [allMax,whichMax] = max(allPercentage);
        
        %0 will be undefined cells. 1 = CD3, 2 = CD11c, 3 = FXIIIA
        phenotype = (allMax>20).*whichMax;
        
        %Get position in phenotype vector where each cell type is
        cd3Idx = find(phenotype==1);
        cd11cIdx = find(phenotype==2);
        fxiiiaIdx = find(phenotype==3);
        
        %Checks labelled mask values for matching values in CD3Idx
        %e.g. cell 25 has value 25 and is kept in mask if CD3Idx contains
        %the value 25. Essentially a matching function
        cd3NucMask = ismember(nuclearMaskLab, cd3Idx);
        cd11cNucMask = ismember(nuclearMaskLab, cd11cIdx);
        fxiiiaNucMask = ismember(nuclearMaskLab, fxiiiaIdx);
        otherNucMask = nuclearMask .* ~cd3NucMask .* ~cd11cNucMask .* ~fxiiiaNucMask;
        
        imwrite(uint8(cd3NucMask*255), [workDir filesep directoryNames{i} filesep 'cd3NucMask.tif']);
        imwrite(uint8(cd11cNucMask*255), [workDir filesep directoryNames{i} filesep 'cd11cNucMask.tif']);
        imwrite(uint8(fxiiiaNucMask*255), [workDir filesep directoryNames{i} filesep 'fxiiiaNucMask.tif']);
        imwrite(uint8(otherNucMask*255), [workDir filesep directoryNames{i} filesep 'otherNucMask.tif']);
        
        %% Perform watershed
        %Use Nick's fancy expandNucleus function to expand out from nuclei
        %to fill cell body. Outputs labelled mask.
        cd3CellMask = expandNucleus(cd3NucMask,cd3Mask);
        cd11cCellMask = expandNucleus(cd11cNucMask,cd11cMask);
        fxiiiaCellMask = expandNucleus(fxiiiaNucMask,fxiiiaMask);
        
        imwrite(uint16(cd3CellMask), [workDir filesep directoryNames{i} filesep 'cd3CellMask.tif']);
        imwrite(uint16(fxiiiaCellMask), [workDir filesep directoryNames{i} filesep 'fxiiiaCellMask.tif']);
        imwrite(uint16(cd11cCellMask), [workDir filesep directoryNames{i} filesep 'cd11cCellMask.tif']);
        imwrite(uint16(otherNucMask), [workDir filesep directoryNames{i} filesep 'otherCellMask.tif']);
%     catch e
%         disp(e.identifier);
%         disp(e.message);
%     end
end