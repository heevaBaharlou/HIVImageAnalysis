%This script will generate csv files for analysis
%Inputs: expandedCellMasks, nucleiMasks, all channels, AF masks, filter
%masks, compMasks, thresholds.csv, membrane Masks, HIV spots mask
%(unfiltered)
%Dependencies: csvwrite_with_thresholds
%Outputs: csvs, HIV spots filtered mask, CD4 Mask


workDir = pwd;
files = dir;
directoryNames = {files([files.isdir]).name};
directoryNames = directoryNames(~ismember(directoryNames,{'.','..'}));

%Must contain thresholds for CD4 and HIV. Other markers have already been
%thresholded and masks created to define cell types. 
[num,txt,raw] = xlsread('thresholds.csv');

[r,~] = size(txt);
fileNames = cell(1,r-1);

for j = 2:r
    fileNames{j-1} = txt{j,1};
end

toRun = [1];

for i = toRun %1:length(directoryNames)
    disp(i)
    disp(directoryNames(i))
%     try
        %% Cell and nuclear masks
        % Nuclear masks
        cd3NucMask = logical(imread([workDir filesep directoryNames{i} filesep 'cd3NucMask.tif']));
        cd3NucMask = double(bwpropfilt(cd3NucMask,'EquivDiameter',[10 50]));
        cd11cNucMask = logical(imread([workDir filesep directoryNames{i} filesep 'cd11cNucMask.tif']));
        cd11cNucMask = bwpropfilt(cd11cNucMask,'EquivDiameter',[10 50]);
        fxiiiaNucMask = logical(imread([workDir filesep directoryNames{i} filesep 'fxiiiaNucMask.tif']));
        fxiiiaNucMask = bwpropfilt(fxiiiaNucMask,'EquivDiameter',[10 50]);
        
        % Cell masks
        cd3CellMask = double(imread([workDir filesep directoryNames{i} filesep 'cd3CellMask.tif']));
        cd11cCellMask = double(imread([workDir filesep directoryNames{i} filesep 'cd11cCellMask.tif']));
        fxiiiaCellMask = double(imread([workDir filesep directoryNames{i} filesep 'fxiiiaCellMask.tif']));
        otherCellMask = logical(imread([workDir filesep directoryNames{i} filesep 'otherNucMask.tif']));
        otherCellMask = bwlabel(bwpropfilt(otherCellMask,'EquivDiameter',[10 50]));
        
        % Label nuclear masks
        %Because expanded segmentation can alter cell number ordering,
        %multiply by 'CellMask' to get Nuclei ordering in same order as
        %'CellMask'
        cd3NucMask = cd3CellMask.*(cd3NucMask>0);
        cd11cNucMask = cd11cCellMask.*(cd11cNucMask>0);
        fxiiiaNucMask = fxiiiaCellMask.*(fxiiiaNucMask>0);
        
        %% Morphology measurements
        % Nuclear Area
        cd3NucArea = regionprops(cd3NucMask, 'Area');
        cd3NucArea = cell2mat(struct2cell(cd3NucArea));
        cd11cNucArea = regionprops(cd11cNucMask, 'Area');
        cd11cNucArea = cell2mat(struct2cell(cd11cNucArea));
        fxiiiaNucArea = regionprops(fxiiiaNucMask, 'Area');
        fxiiiaNucArea = cell2mat(struct2cell(fxiiiaNucArea));
        
        % Cell Area
        cd3CellArea = regionprops(cd3CellMask, 'Area');
        cd3CellArea = cell2mat(struct2cell(cd3CellArea));
        cd11cCellArea = regionprops(cd11cCellMask, 'Area');
        cd11cCellArea = cell2mat(struct2cell(cd11cCellArea));
        fxiiiaCellArea = regionprops(fxiiiaCellMask, 'Area');
        fxiiiaCellArea = cell2mat(struct2cell(fxiiiaCellArea));
        otherCellArea = regionprops(otherCellMask, 'Area');
        otherCellArea = cell2mat(struct2cell(otherCellArea));
        
        % Centroids
       %Always use this weird 2 step way to get centroids
        cd3Centroid = regionprops(cd3NucMask, 'Centroid');
        cd3Centroid = cat(1, cd3Centroid.Centroid);
        cd11cCentroid = regionprops(cd11cNucMask, 'Centroid');
        cd11cCentroid = cat(1, cd11cCentroid.Centroid);
        fxiiiaCentroid = regionprops(fxiiiaNucMask, 'Centroid');
        fxiiiaCentroid = cat(1, fxiiiaCentroid.Centroid);
        otherCentroid = regionprops(otherCellMask, 'Centroid');
        otherCentroid = cat(1, otherCentroid.Centroid);
        
        % Final vector
        nucArea = [cd3NucArea cd11cNucArea fxiiiaNucArea otherCellArea];
        cellArea = [cd3CellArea cd11cCellArea fxiiiaCellArea otherCellArea];
        centroid = [cd3Centroid; cd11cCentroid; fxiiiaCentroid; otherCentroid];
        
        %% Marker intensity
        % Read in markers
        cd3 = double(imread([workDir filesep directoryNames{i} filesep 'CD3_R.tif']));
        cd4 = double(imread([workDir filesep directoryNames{i} filesep 'CD4.tif']));
        cd11c = double(imread([workDir filesep directoryNames{i} filesep 'CD11c.tif']));
        fxiiia = double(imread([workDir filesep directoryNames{i} filesep 'FXIIIA.tif']));
        
        % Read in masks
        cd3Mask = logical(imread([workDir filesep directoryNames{i} filesep 'cd3Mask.tif']));
        cd11cMask = logical(imread([workDir filesep directoryNames{i} filesep 'cd11cMask.tif']));
        fxiiiaMask = logical(imread([workDir filesep directoryNames{i} filesep 'fxiiiaMask.tif']));
        
        %Make and save CD4 Mask
        cd4Mask = imgaussfilt(cd4,1.5);
        
        cd4Thresh = 0;
        
        for j = 1:length(fileNames)
            if strcmp(fileNames{j},directoryNames{i})
                cd4Thresh = num(j,6);
            end
        end
        
        cd4Mask = cd4Mask > cd4Thresh;
        
        imwrite(uint8(cd4Mask)*255, [workDir filesep directoryNames{i} filesep 'cd4Mask.tif']);
        
        % Remove AF and filtered areas
        cd4AF = logical(imread([workDir filesep directoryNames{i} filesep 'CD4_AFidE.tif']));
        cd11cAF = logical(imread([workDir filesep directoryNames{i} filesep 'CD11c_AFidE.tif']));
        fxiiiaAF = logical(imread([workDir filesep directoryNames{i} filesep 'FXIIIA_AFidE.tif']));
        
        %Assign pixels in cd4 image to 0 if that entry is valued 1 in cd4AF image
        cd4(cd4AF) = 0;
        cd4Mask(cd4AF) = 0;
        cd11c(cd11cAF) = 0;
        fxiiia(fxiiiaAF) = 0;
        
        %Can use isFile function instead of try catch
        try
            cd3Filter = logical(imread([workDir filesep directoryNames{i} filesep 'CD3_R_filter.tif']));
            cd3(cd11cFilter) = 0;
        catch
            disp('No CD3 Filter')
        end
        
        try
            cd4Filter = logical(imread([workDir filesep directoryNames{i} filesep 'CD4_filter.tif']));
            cd4(cd4Filter) = 0;
            cd4Mask(cd4Filter) = 0;
        catch
            disp('No CD4 Filter')
        end
        
        try
            cd11cFilter = logical(imread([workDir filesep directoryNames{i} filesep 'CD11c_filter.tif']));
            cd11c(cd11cFilter) = 0;
        catch
            disp('No CD11c Filter')
        end
        
        %Set to NaN so not included in 'mean' measurements
        cd3(cd3==0) = NaN;
        cd4(cd4==0) = NaN;
        cd11c(cd11c==0) = NaN;
        fxiiia(fxiiia==0) = NaN;
        
        % Get nuclear overlap with cell mask
        cd3Nuc_cd3Val = cell2mat(struct2cell(regionprops(cd3NucMask, cd3Mask, 'MeanIntensity'))).*100;
        cd11cNuc_cd3Val = cell2mat(struct2cell(regionprops(cd11cNucMask, cd3Mask, 'MeanIntensity'))).*100;
        fxiiiaNuc_cd3Val = cell2mat(struct2cell(regionprops(fxiiiaNucMask, cd3Mask, 'MeanIntensity'))).*100;
        otherNuc_cd3Val = cell2mat(struct2cell(regionprops(otherCellMask, cd3Mask, 'MeanIntensity'))).*100;
        
        cd3Nuc_cd4Val = cell2mat(struct2cell(regionprops(cd3NucMask, cd4Mask, 'MeanIntensity'))).*100;
        cd11cNuc_cd4Val = cell2mat(struct2cell(regionprops(cd11cNucMask, cd4Mask, 'MeanIntensity'))).*100;
        fxiiiaNuc_cd4Val = cell2mat(struct2cell(regionprops(fxiiiaNucMask, cd4Mask, 'MeanIntensity'))).*100;
        otherNuc_cd4Val = cell2mat(struct2cell(regionprops(otherCellMask, cd4Mask, 'MeanIntensity'))).*100;
        
        cd3Nuc_cd11cVal = cell2mat(struct2cell(regionprops(cd3NucMask, cd11cMask, 'MeanIntensity'))).*100;
        cd11cNuc_cd11cVal = cell2mat(struct2cell(regionprops(cd11cNucMask, cd11cMask, 'MeanIntensity'))).*100;
        fxiiiaNuc_cd11cVal = cell2mat(struct2cell(regionprops(fxiiiaNucMask, cd11cMask, 'MeanIntensity'))).*100;
        otherNuc_cd11cVal = cell2mat(struct2cell(regionprops(otherCellMask, cd11cMask, 'MeanIntensity'))).*100;
        
        cd3Nuc_fxiiiaVal = cell2mat(struct2cell(regionprops(cd3NucMask, fxiiiaMask, 'MeanIntensity'))).*100;
        cd11cNuc_fxiiiaVal = cell2mat(struct2cell(regionprops(cd11cNucMask, fxiiiaMask, 'MeanIntensity'))).*100;
        fxiiiaNuc_fxiiiaVal = cell2mat(struct2cell(regionprops(fxiiiaNucMask, fxiiiaMask, 'MeanIntensity'))).*100;
        otherNuc_fxiiiaVal = cell2mat(struct2cell(regionprops(otherCellMask, fxiiiaMask, 'MeanIntensity'))).*100;
        
        % Get cell mean intensity
        %PixelValues gets list of pixel values for each cell
        %After converting to 'cell' then use 'lapply' equivalent 'cellfun'
        %to apply nanmean to each cell. '@' used as function prefix in
        %cellfun.
        cd3Cell_cd3Val = cellfun(@nanmean, struct2cell(regionprops(cd3CellMask, cd3, 'PixelValues')));
        %If whole cell is NaN's then set to 0.
        cd3Cell_cd3Val(isnan(cd3Cell_cd3Val)) = 0;
        cd11cCell_cd3Val = cellfun(@nanmean, struct2cell(regionprops(cd11cCellMask, cd3, 'PixelValues')));
        cd11cCell_cd3Val(isnan(cd11cCell_cd3Val)) = 0;
        fxiiiaCell_cd3Val = cellfun(@nanmean, struct2cell(regionprops(fxiiiaCellMask, cd3, 'PixelValues')));
        fxiiiaCell_cd3Val(isnan(fxiiiaCell_cd3Val)) = 0;
        otherCell_cd3Val = cellfun(@nanmean, struct2cell(regionprops(otherCellMask, cd3, 'PixelValues')));
        otherCell_cd3Val(isnan(otherCell_cd3Val)) = 0;
        
        cd3Cell_cd4Val = cellfun(@nanmean, struct2cell(regionprops(cd3CellMask, cd4, 'PixelValues')));
        cd3Cell_cd4Val(isnan(cd3Cell_cd4Val)) = 0;
        cd11cCell_cd4Val = cellfun(@nanmean, struct2cell(regionprops(cd11cCellMask, cd4, 'PixelValues')));
        cd11cCell_cd4Val(isnan(cd11cCell_cd4Val)) = 0;
        fxiiiaCell_cd4Val = cellfun(@nanmean, struct2cell(regionprops(fxiiiaCellMask, cd4, 'PixelValues')));
        fxiiiaCell_cd4Val(isnan(fxiiiaCell_cd4Val)) = 0;
        otherCell_cd4Val = cellfun(@nanmean, struct2cell(regionprops(otherCellMask, cd4, 'PixelValues')));
        otherCell_cd4Val(isnan(otherCell_cd4Val)) = 0;
        
        cd3Cell_cd11cVal = cellfun(@nanmean, struct2cell(regionprops(cd3CellMask, cd11c, 'PixelValues')));
        cd3Cell_cd11cVal(isnan(cd3Cell_cd11cVal)) = 0;
        cd11cCell_cd11cVal = cellfun(@nanmean, struct2cell(regionprops(cd11cCellMask, cd11c, 'PixelValues')));
        cd11cCell_cd11cVal(isnan(cd11cCell_cd11cVal)) = 0;
        fxiiiaCell_cd11cVal = cellfun(@nanmean, struct2cell(regionprops(fxiiiaCellMask, cd11c, 'PixelValues')));
        fxiiiaCell_cd11cVal(isnan(fxiiiaCell_cd11cVal)) = 0;
        otherCell_cd11cVal = cellfun(@nanmean, struct2cell(regionprops(otherCellMask, cd11c, 'PixelValues')));
        otherCell_cd11cVal(isnan(otherCell_cd11cVal)) = 0;
        
        cd3Cell_fxiiiaVal = cellfun(@nanmean, struct2cell(regionprops(cd3CellMask, fxiiia, 'PixelValues')));
        cd3Cell_fxiiiaVal(isnan(cd3Cell_fxiiiaVal)) = 0;
        cd11cCell_fxiiiaVal = cellfun(@nanmean, struct2cell(regionprops(cd11cCellMask, fxiiia, 'PixelValues')));
        cd11cCell_fxiiiaVal(isnan(cd11cCell_fxiiiaVal)) = 0;
        fxiiiaCell_fxiiiaVal = cellfun(@nanmean, struct2cell(regionprops(fxiiiaCellMask, fxiiia, 'PixelValues')));
        fxiiiaCell_fxiiiaVal(isnan(fxiiiaCell_fxiiiaVal)) = 0;
        otherCell_fxiiiaVal = cellfun(@nanmean, struct2cell(regionprops(otherCellMask, fxiiia, 'PixelValues')));
        otherCell_fxiiiaVal(isnan(otherCell_fxiiiaVal)) = 0;
        
        % Final vector
        cd3Nuc = [cd3Nuc_cd3Val cd11cNuc_cd3Val fxiiiaNuc_cd3Val otherNuc_cd3Val];
        cd3Cell = [cd3Cell_cd3Val cd11cCell_cd3Val fxiiiaCell_cd3Val otherCell_cd3Val];
        
        cd4Nuc = [cd3Nuc_cd4Val cd11cNuc_cd4Val fxiiiaNuc_cd4Val otherNuc_cd4Val];
        cd4Cell = [cd3Cell_cd4Val cd11cCell_cd4Val fxiiiaCell_cd4Val otherCell_cd4Val];
        
        cd11cNuc = [cd3Nuc_cd11cVal cd11cNuc_cd11cVal fxiiiaNuc_cd11cVal otherNuc_cd11cVal];
        cd11cCell = [cd3Cell_cd11cVal cd11cCell_cd11cVal fxiiiaCell_cd11cVal otherCell_cd11cVal];
        
        fxiiiaNuc = [cd3Nuc_fxiiiaVal cd11cNuc_fxiiiaVal fxiiiaNuc_fxiiiaVal otherNuc_fxiiiaVal];
        fxiiiaCell = [cd3Cell_fxiiiaVal cd11cCell_fxiiiaVal fxiiiaCell_fxiiiaVal otherCell_fxiiiaVal];
        
        nCell = length(fxiiiaCell);
        
        %% Free up memory
        clearvars cd3 cd4 cd11c fxiiia
        clearvars cd3Mask cd4Mask cd11cMask fxiiiaMask 
        
        %% HIV
        %Filter previously generated hiv spots mask using thresholds and
        %removing AF. Then get HIV area, spot number and distance for each
        %cell
        hivArea = zeros(1,nCell);
        hivSpotCount = zeros(1,nCell);
        hivDist = zeros(1,nCell);
        
        if contains(lower(directoryNames{i}),'hiv')
            % Read in images
            hiv = double(imread([workDir filesep directoryNames{i} filesep 'HIVRNA.tif']));
            hivSpots = logical(imread([workDir filesep directoryNames{i} filesep 'spots.tif']));
            hivAF = logical(imread([workDir filesep directoryNames{i} filesep 'HIVRNA_AFidE.tif']));
            
            hivMask = imgaussfilt(hiv,1.5);
            
            hivThresh = 0;
            
            for j = 1:length(fileNames)
                if strcmp(fileNames{j},directoryNames{i})
                    hivThresh = num(j,5);
                end
            end
            
            hivMask = hivMask > hivThresh;
            hivMask = (hivMask - hivAF)>0;
            
            hivSpots = hivSpots .* hivMask;
            
            imwrite(uint8(hivMask)*255, [workDir filesep directoryNames{i} filesep 'hivMask.tif']);
            imwrite(uint8(hivSpots)*255, [workDir filesep directoryNames{i} filesep 'spotsFiltered.tif']);
            
            % HIV Area
            cd3_hivVal = cell2mat(struct2cell(regionprops(cd3CellMask, hivMask, 'MeanIntensity'))).*cd3CellArea;
            cd11c_hivVal = cell2mat(struct2cell(regionprops(cd11cCellMask, hivMask, 'MeanIntensity'))).*cd11cCellArea;
            fxiiia_hivVal = cell2mat(struct2cell(regionprops(fxiiiaCellMask, hivMask, 'MeanIntensity'))).*fxiiiaCellArea;
            other_hivVal = cell2mat(struct2cell(regionprops(otherCellMask, hivMask, 'MeanIntensity'))).*otherCellArea;
            
            % HIV Spots
            cd3_hivCount = cell2mat(struct2cell(regionprops(cd3CellMask, hivSpots, 'MeanIntensity'))).*cd3CellArea;
            cd11c_hivCount = cell2mat(struct2cell(regionprops(cd11cCellMask, hivSpots, 'MeanIntensity'))).*cd11cCellArea;
            fxiiia_hivCount = cell2mat(struct2cell(regionprops(fxiiiaCellMask, hivSpots, 'MeanIntensity'))).*fxiiiaCellArea;
            other_hivCount = cell2mat(struct2cell(regionprops(otherCellMask, hivSpots, 'MeanIntensity'))).*otherCellArea;
            
            % Distance map
            hivDistMap = bwdist(hivSpots);
            cd3_hivDist = cell2mat(struct2cell(regionprops(cd3CellMask, hivDistMap, 'MinIntensity')));
            cd11c_hivDist = cell2mat(struct2cell(regionprops(cd11cCellMask, hivDistMap, 'MinIntensity')));
            fxiiia_hivDist = cell2mat(struct2cell(regionprops(fxiiiaCellMask, hivDistMap, 'MinIntensity')));
            other_hivDist = cell2mat(struct2cell(regionprops(otherCellMask, hivDistMap, 'MinIntensity')));
            
            % Final Vectors
            hivArea = [cd3_hivVal cd11c_hivVal fxiiia_hivVal other_hivVal];
            hivSpotCount = [cd3_hivCount cd11c_hivCount fxiiia_hivCount other_hivCount];
            hivDist = [cd3_hivDist cd11c_hivDist fxiiia_hivDist other_hivDist];
        end
        
  
        
        %% Free up memory
        clearvars cd3CellArea cd11cCellArea fxiiiaCellArea otherCellArea
        clearvars hiv hivAF hivMask  cd4AF cd11cAF fxiiiaAF cd11cNucMask cd3NucMask fxiiiaNucMask
        clearvars cd3_hivVal cd11c_hivVal fxiiia_hivVal other_hivVal
        clearvars cd3_hivCount cd11c_hivCount fxiiia_hivCount other_hivCount
        clearvars cd3_hivDist cd11c_hivDist fxiiia_hivDist other_hivDist
        
        %% Compartments
        
        %pre-allocation of vectors which will become distance of each cell
        %from each compartment
        epDist = zeros(1,nCell);
        laDist = zeros(1,nCell);
        lpDist = zeros(1,nCell);
        smDist = zeros(1,nCell);
        laInnerDist = zeros(1,nCell);
        
        lpComp = zeros(1,nCell);
        smComp = zeros(1,nCell);
        epComp = zeros(1,nCell);
        laComp = zeros(1,nCell);
        %Defined as ones. Changed to zero later if it is in a defined
        %compartment. In loop below.
        unComp = ones(1,nCell);
        
        lpMask = zeros(size(otherCellMask));
        smMask = zeros(size(otherCellMask));
        epMask = zeros(size(otherCellMask));
        laMask = zeros(size(otherCellMask));
        
        epDistMask = zeros(size(otherCellMask));
        laDistMask = zeros(size(otherCellMask));
        lpDistMask = zeros(size(otherCellMask));
        
        compMask = imread([workDir filesep directoryNames{i} filesep 'compMask.tif']);
        lpMask(compMask == 4) = 1;
        smMask(compMask == 2) = 1;
        epMask(compMask == 1) = 1;
        laMask(compMask == 3) = 1;
        
        %Measure percentage of each cellType in each compartment
        
        cd3_lpVal = cell2mat(struct2cell(regionprops(cd3CellMask, lpMask, 'MeanIntensity')));
        cd11c_lpVal = cell2mat(struct2cell(regionprops(cd11cCellMask, lpMask, 'MeanIntensity')));
        fxiiia_lpVal = cell2mat(struct2cell(regionprops(fxiiiaCellMask, lpMask, 'MeanIntensity')));
        other_lpVal = cell2mat(struct2cell(regionprops(otherCellMask, lpMask, 'MeanIntensity')));
        
        lpVal = [cd3_lpVal cd11c_lpVal fxiiia_lpVal other_lpVal];
        clearvars cd3_lpVal cd11c_lpVal fxiiia_lpVal other_lpVal
        
        cd3_smVal = cell2mat(struct2cell(regionprops(cd3CellMask, smMask, 'MeanIntensity')));
        cd11c_smVal = cell2mat(struct2cell(regionprops(cd11cCellMask, smMask, 'MeanIntensity')));
        fxiiia_smVal = cell2mat(struct2cell(regionprops(fxiiiaCellMask, smMask, 'MeanIntensity')));
        other_smVal = cell2mat(struct2cell(regionprops(otherCellMask, smMask, 'MeanIntensity')));
        
        smVal = [cd3_smVal cd11c_smVal fxiiia_smVal other_smVal];
        clearvars cd3_smVal cd11c_smVal fxiiia_smVal other_smVal
        
        cd3_epVal = cell2mat(struct2cell(regionprops(cd3CellMask, epMask, 'MeanIntensity')));
        cd11c_epVal = cell2mat(struct2cell(regionprops(cd11cCellMask, epMask, 'MeanIntensity')));
        fxiiia_epVal = cell2mat(struct2cell(regionprops(fxiiiaCellMask, epMask, 'MeanIntensity')));
        other_epVal = cell2mat(struct2cell(regionprops(otherCellMask, epMask, 'MeanIntensity')));
        
        epVal = [cd3_epVal cd11c_epVal fxiiia_epVal other_epVal];
        clearvars cd3_epVal cd11c_epVal fxiiia_epVal other_epVal
        
        cd3_laVal = cell2mat(struct2cell(regionprops(cd3CellMask, laMask, 'MeanIntensity')));
        cd11c_laVal = cell2mat(struct2cell(regionprops(cd11cCellMask, laMask, 'MeanIntensity')));
        fxiiia_laVal = cell2mat(struct2cell(regionprops(fxiiiaCellMask, laMask, 'MeanIntensity')));
        other_laVal = cell2mat(struct2cell(regionprops(otherCellMask, laMask, 'MeanIntensity')));
        
        laVal = [cd3_laVal cd11c_laVal fxiiia_laVal other_laVal];
        clearvars cd3_epVal cd11c_epVal fxiiia_epVal other_epVal
        
        %Assign cells to compartments in a hierarchy with LAs having highest priority.  
        for j = 1:nCell
            if laVal(j) > 0.25
                laComp(j) = 1;  
                unComp(j) = 0;
            elseif smVal(j) > 0.25
                smComp(j) = 1;
                unComp(j) = 0;
            elseif epVal(j) > 0.25
                epComp(j) = 1;
                unComp(j) = 0;
            elseif lpVal(j) > 0.25
                lpComp(j) = 1;
                unComp(j) = 0;
            end
        end
        
        %If compartment exists then create distance map and measure
        %distance for each cell in each cell type. 
        if max(epMask(:)) > 0
            epDistMask = bwdist(epMask);
            
            cd3_epDist = cell2mat(struct2cell(regionprops(cd3CellMask, epDistMask, 'MeanIntensity')));
            cd11c_epDist = cell2mat(struct2cell(regionprops(cd11cCellMask, epDistMask, 'MeanIntensity')));
            fxiiia_epDist = cell2mat(struct2cell(regionprops(fxiiiaCellMask, epDistMask, 'MeanIntensity')));
            other_epDist = cell2mat(struct2cell(regionprops(otherCellMask, epDistMask, 'MeanIntensity')));
            
            epDist = [cd3_epDist cd11c_epDist fxiiia_epDist other_epDist];
             clearvars epDistMask cd3_epDist cd11c_epDist fxiiia_epDist other_epDist
        end
        
        if max(laMask(:)) > 0
            laDistMask = bwdist(laMask);
            
            cd3_laDist = cell2mat(struct2cell(regionprops(cd3CellMask, laDistMask, 'MeanIntensity')));
            cd11c_laDist = cell2mat(struct2cell(regionprops(cd11cCellMask, laDistMask, 'MeanIntensity')));
            fxiiia_laDist = cell2mat(struct2cell(regionprops(fxiiiaCellMask, laDistMask, 'MeanIntensity')));
            other_laDist = cell2mat(struct2cell(regionprops(otherCellMask, laDistMask, 'MeanIntensity')));
            
            laDist = [cd3_laDist cd11c_laDist fxiiia_laDist other_laDist];
             clearvars laDistMask cd3_laDist cd11c_laDist fxiiia_laDist other_laDist
        end
        
        
        if max(lpMask(:)) > 0
            lpDistMask = bwdist(lpMask);
            
            cd3_lpDist = cell2mat(struct2cell(regionprops(cd3CellMask, lpDistMask, 'MeanIntensity')));
            cd11c_lpDist = cell2mat(struct2cell(regionprops(cd11cCellMask, lpDistMask, 'MeanIntensity')));
            fxiiia_lpDist = cell2mat(struct2cell(regionprops(fxiiiaCellMask, lpDistMask, 'MeanIntensity')));
            other_lpDist = cell2mat(struct2cell(regionprops(otherCellMask, lpDistMask, 'MeanIntensity')));
            
            lpDist = [cd3_lpDist cd11c_lpDist fxiiia_lpDist other_lpDist];
             clearvars lpDistMask cd3_lpDist cd11c_lpDist fxiiia_lpDist other_lpDist
        end
        
        
         if max(smMask(:)) > 0
            smDistMask = bwdist(smMask);
            
            cd3_smDist = cell2mat(struct2cell(regionprops(cd3CellMask, smDistMask, 'MeanIntensity')));
            cd11c_smDist = cell2mat(struct2cell(regionprops(cd11cCellMask, smDistMask, 'MeanIntensity')));
            fxiiia_smDist = cell2mat(struct2cell(regionprops(fxiiiaCellMask, smDistMask, 'MeanIntensity')));
            other_smDist = cell2mat(struct2cell(regionprops(otherCellMask, smDistMask, 'MeanIntensity')));
            
            smDist = [cd3_smDist cd11c_smDist fxiiia_smDist other_smDist];
            clearvars smDistMask cd3_smDist cd11c_smDist fxiiia_smDist other_smDist
        end
        
        
        % Get LAInnerDist (Distance map from LA edge emanating inward).
        if max(laMask(:)) > 0
            laMask(compMask ~= 3) = 1;
            laMask(compMask == 3) = 0;
            laDistMask = bwdist(laMask);
            
            cd3_laDist = cell2mat(struct2cell(regionprops(cd3CellMask, laDistMask, 'MeanIntensity')));
            cd11c_laDist = cell2mat(struct2cell(regionprops(cd11cCellMask, laDistMask, 'MeanIntensity')));
            fxiiia_laDist = cell2mat(struct2cell(regionprops(fxiiiaCellMask, laDistMask, 'MeanIntensity')));
            other_laDist = cell2mat(struct2cell(regionprops(otherCellMask, laDistMask, 'MeanIntensity')));
            
            laInnerDist = [cd3_laDist cd11c_laDist fxiiia_laDist other_laDist];
            clearvars laDistMask cd3_laDist cd11c_laDist fxiiia_laDist other_laDist
        end
        
                
        %% Compartment HIV distances
        
        %Variables that will hold cell distances from HIV in a compartment.
        %Eg lpHIVDist would be the nearest distance of any cell in the
        %image to LP.
        
        lpHIVDist = zeros(1,nCell);
        smHIVDist = zeros(1,nCell);
        epHIVDist = zeros(1,nCell);
        laHIVDist = zeros(1,nCell);       
        
        % This part was already defined in above section on compartments
%         lpMask = zeros(size(otherCellMask));
%         smMask = zeros(size(otherCellMask));
%         epMask = zeros(size(otherCellMask));
%         laMask = zeros(size(otherCellMask));
% 
%         compMask = imread([workDir filesep directoryNames{i} filesep 'compMask.tif']);
%         lpMask(compMask == 4) = 1;
%         smMask(compMask == 2) = 1;
%         epMask(compMask == 1) = 1;
%         laMask(compMask == 3) = 1;
        
        %Create masks of HIV for each compartment
        lpHIVSpots = lpMask & hivSpots;
        smHIVSpots = smMask & hivSpots;
        epHIVSpots = epMask & hivSpots;
        laHIVSpots = laMask & hivSpots;
        
        %Get dist of cells from HIV+ cells from each compartment
        
        % LP
        if max(lpMask(:)) > 0
            hivDistMap = bwdist(lpHIVSpots);
            cd3_hivDist = cell2mat(struct2cell(regionprops(cd3CellMask, hivDistMap, 'MinIntensity')));
            cd11c_hivDist = cell2mat(struct2cell(regionprops(cd11cCellMask, hivDistMap, 'MinIntensity')));
            fxiiia_hivDist = cell2mat(struct2cell(regionprops(fxiiiaCellMask, hivDistMap, 'MinIntensity')));
            other_hivDist = cell2mat(struct2cell(regionprops(otherCellMask, hivDistMap, 'MinIntensity')));
            %final vector
            lpHIVDist = [cd3_hivDist cd11c_hivDist fxiiia_hivDist other_hivDist];
            clearvars cd3_hivDist cd11c_hivDist fxiiia_hivDist other_hivDist
        end
        
        % SM
        if max(smMask(:)) > 0
            hivDistMap = bwdist(smHIVSpots);
            cd3_hivDist = cell2mat(struct2cell(regionprops(cd3CellMask, hivDistMap, 'MinIntensity')));
            cd11c_hivDist = cell2mat(struct2cell(regionprops(cd11cCellMask, hivDistMap, 'MinIntensity')));
            fxiiia_hivDist = cell2mat(struct2cell(regionprops(fxiiiaCellMask, hivDistMap, 'MinIntensity')));
            other_hivDist = cell2mat(struct2cell(regionprops(otherCellMask, hivDistMap, 'MinIntensity')));
            %final vector
            smHIVDist = [cd3_hivDist cd11c_hivDist fxiiia_hivDist other_hivDist];
            clearvars cd3_hivDist cd11c_hivDist fxiiia_hivDist other_hivDist
        end
        

        % EP
        if max(epMask(:)) > 0
            hivDistMap = bwdist(epHIVSpots);
            cd3_hivDist = cell2mat(struct2cell(regionprops(cd3CellMask, hivDistMap, 'MinIntensity')));
            cd11c_hivDist = cell2mat(struct2cell(regionprops(cd11cCellMask, hivDistMap, 'MinIntensity')));
            fxiiia_hivDist = cell2mat(struct2cell(regionprops(fxiiiaCellMask, hivDistMap, 'MinIntensity')));
            other_hivDist = cell2mat(struct2cell(regionprops(otherCellMask, hivDistMap, 'MinIntensity')));
            %final vector
            epHIVDist = [cd3_hivDist cd11c_hivDist fxiiia_hivDist other_hivDist];
            clearvars cd3_hivDist cd11c_hivDist fxiiia_hivDist other_hivDist
        end
        
        %LA
        if max(laMask(:)) > 0
            hivDistMap = bwdist(laHIVSpots);
            cd3_hivDist = cell2mat(struct2cell(regionprops(cd3CellMask, hivDistMap, 'MinIntensity')));
            cd11c_hivDist = cell2mat(struct2cell(regionprops(cd11cCellMask, hivDistMap, 'MinIntensity')));
            fxiiia_hivDist = cell2mat(struct2cell(regionprops(fxiiiaCellMask, hivDistMap, 'MinIntensity')));
            other_hivDist = cell2mat(struct2cell(regionprops(otherCellMask, hivDistMap, 'MinIntensity')));
            %final vector
            laHIVDist = [cd3_hivDist cd11c_hivDist fxiiia_hivDist other_hivDist];
            clearvars cd3_hivDist cd11c_hivDist fxiiia_hivDist other_hivDist
        end

        %% Write csv
        cd3Count = max(cd3CellMask(:));
        cd11cCount = max(cd11cCellMask(:));
        fxiiiaCount = max(fxiiiaCellMask(:));
        otherCount = max(otherCellMask(:));
        
        cellIdx = 1:nCell;
                
        %Stack cellTypes in the spreadsheet with CD3 first.
        cellType = zeros(1,nCell);
        cellType(1:cd3Count) = 1;
        cellType((cd3Count+1):(cd3Count+cd11cCount)) = 2;
        cellType((cd3Count+cd11cCount+1):(cd3Count+cd11cCount+fxiiiaCount)) = 3;
      
        headers = {'CellID', 'CellType', 'X', 'Y', 'NucArea', 'CellArea', ...
            'CD3Overlap', 'CD3', 'CD4Overlap', 'CD4', ...
            'CD11cOverlap', 'CD11c', 'FXIIIaOverlap', 'FXIIIA', ...
            'hivArea', 'hivCount', 'hivDist', ...
            'EP', 'LP', 'LA', 'SM', 'Undefined', ...
            'EPDist', 'LPDist', 'LADist', 'SMDist', 'LAInnerDist', 'lpHIVDist', 'smHIVDist', 'epHIVDist', 'laHIVDist'};
        %Transpose everything except for centroid which is already two
        %column vectors
        data = [cellIdx' cellType' centroid nucArea' cellArea' ...
            cd3Nuc' cd3Cell' cd4Nuc' cd4Cell' ...
            cd11cNuc' cd11cCell' fxiiiaNuc' fxiiiaCell' ...
            hivArea' hivSpotCount' hivDist' ...
            epComp' lpComp' laComp' smComp' unComp' ...
            epDist' lpDist' laDist' smDist' laInnerDist' smHIVDist' smHIVDist' epHIVDist' laHIVDist'];
        
        csvwrite_with_headers([workDir filesep directoryNames{i} filesep 'cellData.csv'], data, headers);
        csvwrite_with_headers([workDir filesep 'cellData' filesep directoryNames{i} '.csv'], data, headers);
        
        %Free up memory
        clearvars -except workDir files directoryNames fileNames toRun num txt raw r i j
        
        

%         %% Neighbours
%         expansion = 18;
%         se = strel('disk', expansion);
%         
%         neighbrs = zeros(nCell,100);
%         currentMax = 0;
%         
%         % CD3
%         [cd3CellId, cd3Cell_cd3Patch, cd3Cell_cd11cPatch, cd3Cell_fxiiiaPatch, cd3Cell_otherPatch] = ...
%             getPatch(cd3CellMask, cd11cCellMask, fxiiiaCellMask, otherCellMask, expansion);
%         cd3Extended = cellfun(@(x,y) imdilate(x==y,se,'same'),cd3Cell_cd3Patch,num2cell(cd3CellId),'UniformOutput',false);
%         
%         cd3_cd3Overlap = cellfun(@(x,y) x(y),cd3Cell_cd3Patch,cd3Extended,'UniformOutput',false);
%         cd3_cd3Neighbrs = cellfun(@(x) unique(x),cd3_cd3Overlap,'UniformOutput',false);
%         cd3_cd11cOverlap = cellfun(@(x,y) x(y),cd3Cell_cd11cPatch,cd3Extended,'UniformOutput',false);
%         cd3_cd11cNeighbrs = cellfun(@(x) unique(x),cd3_cd11cOverlap,'UniformOutput',false);
%         cd3_fxiiiaOverlap = cellfun(@(x,y) x(y),cd3Cell_fxiiiaPatch,cd3Extended,'UniformOutput',false);
%         cd3_fxiiiaNeighbrs = cellfun(@(x) unique(x),cd3_fxiiiaOverlap,'UniformOutput',false);
%         cd3_otherOverlap = cellfun(@(x,y) x(y),cd3Cell_otherPatch,cd3Extended,'UniformOutput',false);
%         cd3_otherNeighbrs = cellfun(@(x) unique(x),cd3_otherOverlap,'UniformOutput',false);
%         
%         for j = 1:cd3Count
%             neighbrs1 = cd3_cd3Neighbrs{j}';
%             neighbrs1 = neighbrs1(neighbrs1~=0 & neighbrs1~=j);
%             
%             neighbrs2 = cd3_cd11cNeighbrs{j}';
%             neighbrs2 = neighbrs2(neighbrs2~=0 & neighbrs2~=j) + cd3Count;
%             
%             neighbrs3 = cd3_fxiiiaNeighbrs{j}';
%             neighbrs3 = neighbrs3(neighbrs3~=0 & neighbrs3~=j) + cd3Count + cd11cCount;
%             
%             neighbrs4 = cd3_otherNeighbrs{j}';
%             neighbrs4 = neighbrs4(neighbrs4~=0 & neighbrs4~=j) + cd3Count + cd11cCount + fxiiiaCount;
%             
%             allNeighbrs = sort([neighbrs1 neighbrs2 neighbrs3 neighbrs4]);
%             
%             if length(allNeighbrs) > currentMax
%                 currentMax = length(allNeighbrs);
%             end
%             
%             neighbrs(j,1:length(allNeighbrs)) = allNeighbrs;
%         end
%         
%         % CD11c
%         [cd11cCellId, cd11cCell_cd11cPatch, cd11cCell_fxiiiaPatch, cd11cCell_otherPatch, cd11cCell_cd3Patch] = ...
%             getPatch(cd11cCellMask, fxiiiaCellMask, otherCellMask, cd3CellMask, expansion);
%         cd11cExtended = cellfun(@(x,y) imdilate(x==y,se,'same'),cd11cCell_cd11cPatch,num2cell(cd11cCellId),'UniformOutput',false);
%         
%         cd11c_cd3Overlap = cellfun(@(x,y) x(y),cd11cCell_cd3Patch,cd11cExtended,'UniformOutput',false);
%         cd11c_cd3Neighbrs = cellfun(@(x) unique(x),cd11c_cd3Overlap,'UniformOutput',false);
%         cd11c_cd11cOverlap = cellfun(@(x,y) x(y),cd11cCell_cd11cPatch,cd11cExtended,'UniformOutput',false);
%         cd11c_cd11cNeighbrs = cellfun(@(x) unique(x),cd11c_cd11cOverlap,'UniformOutput',false);
%         cd11c_fxiiiaOverlap = cellfun(@(x,y) x(y),cd11cCell_fxiiiaPatch,cd11cExtended,'UniformOutput',false);
%         cd11c_fxiiiaNeighbrs = cellfun(@(x) unique(x),cd11c_fxiiiaOverlap,'UniformOutput',false);
%         cd11c_otherOverlap = cellfun(@(x,y) x(y),cd11cCell_otherPatch,cd11cExtended,'UniformOutput',false);
%         cd11c_otherNeighbrs = cellfun(@(x) unique(x),cd11c_otherOverlap,'UniformOutput',false);
%         
%         for j = 1:cd11cCount
%             neighbrs1 = cd11c_cd3Neighbrs{j}';
%             neighbrs1 = neighbrs1(neighbrs1~=0 & neighbrs1~=j);
%             
%             neighbrs2 = cd11c_cd11cNeighbrs{j}';
%             neighbrs2 = neighbrs2(neighbrs2~=0 & neighbrs2~=j) + cd3Count;
%             
%             neighbrs3 = cd11c_fxiiiaNeighbrs{j}';
%             neighbrs3 = neighbrs3(neighbrs3~=0 & neighbrs3~=j) + cd3Count + cd11cCount;
%             
%             neighbrs4 = cd11c_otherNeighbrs{j}';
%             neighbrs4 = neighbrs4(neighbrs4~=0 & neighbrs4~=j) + cd3Count + cd11cCount + fxiiiaCount;
%             
%             allNeighbrs = sort([neighbrs1 neighbrs2 neighbrs3 neighbrs4]);
%             
%             if length(allNeighbrs) > currentMax
%                 currentMax = length(allNeighbrs);
%             end
%             
%             neighbrs(j+cd3Count,1:length(allNeighbrs)) = allNeighbrs;
%         end
%         
%         % FXIIIa
%         [fxiiiaCellId, fxiiiaCell_fxiiiaPatch, fxiiiaCell_otherPatch, fxiiiaCell_cd3Patch, fxiiiaCell_cd11cPatch] = ...
%             getPatch(fxiiiaCellMask, otherCellMask, cd3CellMask, cd11cCellMask, expansion);
%         fxiiiaExtended = cellfun(@(x,y) imdilate(x==y,se,'same'),fxiiiaCell_fxiiiaPatch,num2cell(fxiiiaCellId),'UniformOutput',false);
%         
%         fxiiia_cd3Overlap = cellfun(@(x,y) x(y),fxiiiaCell_cd3Patch,fxiiiaExtended,'UniformOutput',false);
%         fxiiia_cd3Neighbrs = cellfun(@(x) unique(x),fxiiia_cd3Overlap,'UniformOutput',false);
%         fxiiia_cd11cOverlap = cellfun(@(x,y) x(y),fxiiiaCell_cd11cPatch,fxiiiaExtended,'UniformOutput',false);
%         fxiiia_cd11cNeighbrs = cellfun(@(x) unique(x),fxiiia_cd11cOverlap,'UniformOutput',false);
%         fxiiia_fxiiiaOverlap = cellfun(@(x,y) x(y),fxiiiaCell_fxiiiaPatch,fxiiiaExtended,'UniformOutput',false);
%         fxiiia_fxiiiaNeighbrs = cellfun(@(x) unique(x),fxiiia_fxiiiaOverlap,'UniformOutput',false);
%         fxiiia_otherOverlap = cellfun(@(x,y) x(y),fxiiiaCell_otherPatch,fxiiiaExtended,'UniformOutput',false);
%         fxiiia_otherNeighbrs = cellfun(@(x) unique(x),fxiiia_otherOverlap,'UniformOutput',false);
%         
%         for j = 1:fxiiiaCount
%             neighbrs1 = fxiiia_cd3Neighbrs{j}';
%             neighbrs1 = neighbrs1(neighbrs1~=0 & neighbrs1~=j);
%             
%             neighbrs2 = fxiiia_cd11cNeighbrs{j}';
%             neighbrs2 = neighbrs2(neighbrs2~=0 & neighbrs2~=j) + cd3Count;
%             
%             neighbrs3 = fxiiia_fxiiiaNeighbrs{j}';
%             neighbrs3 = neighbrs3(neighbrs3~=0 & neighbrs3~=j) + cd3Count + cd11cCount;
%             
%             neighbrs4 = fxiiia_otherNeighbrs{j}';
%             neighbrs4 = neighbrs4(neighbrs4~=0 & neighbrs4~=j) + cd3Count + cd11cCount + fxiiiaCount;
%             
%             allNeighbrs = sort([neighbrs1 neighbrs2 neighbrs3 neighbrs4]);
%             
%             if length(allNeighbrs) > currentMax
%                 currentMax = length(allNeighbrs);
%             end
%             
%             neighbrs(j+cd3Count+cd11cCount,1:length(allNeighbrs)) = allNeighbrs;
%         end
%         
%         % Other
%         [otherCellId, otherCell_otherPatch, otherCell_cd3Patch, otherCell_cd11cPatch, otherCell_fxiiiaPatch] = ...
%             getPatch(otherCellMask, cd3CellMask, cd11cCellMask, fxiiiaCellMask, expansion);
%         otherExtended = cellfun(@(x,y) imdilate(x==y,se,'same'),otherCell_fxiiiaPatch,num2cell(otherCellId),'UniformOutput',false);
%         
%         other_cd3Overlap = cellfun(@(x,y) x(y),otherCell_cd3Patch,otherExtended,'UniformOutput',false);
%         other_cd3Neighbrs = cellfun(@(x) unique(x),other_cd3Overlap,'UniformOutput',false);
%         other_cd11cOverlap = cellfun(@(x,y) x(y),otherCell_cd11cPatch,otherExtended,'UniformOutput',false);
%         other_cd11cNeighbrs = cellfun(@(x) unique(x),other_cd11cOverlap,'UniformOutput',false);
%         other_fxiiiaOverlap = cellfun(@(x,y) x(y),otherCell_fxiiiaPatch,otherExtended,'UniformOutput',false);
%         other_fxiiiaNeighbrs = cellfun(@(x) unique(x),other_fxiiiaOverlap,'UniformOutput',false);
%         other_otherOverlap = cellfun(@(x,y) x(y),otherCell_otherPatch,otherExtended,'UniformOutput',false);
%         other_otherNeighbrs = cellfun(@(x) unique(x),other_otherOverlap,'UniformOutput',false);
%         
%         for j = 1:fxiiiaCount
%             neighbrs1 = other_cd3Neighbrs{j}';
%             neighbrs1 = neighbrs1(neighbrs1~=0 & neighbrs1~=j);
%             
%             neighbrs2 = other_cd11cNeighbrs{j}';
%             neighbrs2 = neighbrs2(neighbrs2~=0 & neighbrs2~=j) + cd3Count;
%             
%             neighbrs3 = other_fxiiiaNeighbrs{j}';
%             neighbrs3 = neighbrs3(neighbrs3~=0 & neighbrs3~=j) + cd3Count + cd11cCount;
%             
%             neighbrs4 = other_otherNeighbrs{j}';
%             neighbrs4 = neighbrs4(neighbrs4~=0 & neighbrs4~=j) + cd3Count + cd11cCount + fxiiiaCount;
%             
%             allNeighbrs = sort([neighbrs1 neighbrs2 neighbrs3 neighbrs4]);
%             
%             if length(allNeighbrs) > currentMax
%                 currentMax = length(allNeighbrs);
%             end
%             
%             neighbrs(j+cd3Count+cd11cCount+fxiiiaCount,1:length(allNeighbrs)) = allNeighbrs;
%         end
%         
%         neighbrHeaders =  cellstr(num2str((0:currentMax)'))';
%         neighbrHeaders = strrep(neighbrHeaders,' ','');
%         neighbrHeaders = strcat('N', neighbrHeaders);
%         neighbrHeaders{1} = 'CellId';
%         neighbrsData = [cellIdx' neighbrs(:,1:currentMax)];
%         
%         csvwrite_with_headers([workDir filesep directoryNames{i} filesep 'neighbrs.csv'], neighbrsData, neighbrHeaders);
%         csvwrite_with_headers([workDir filesep 'neighbrs' filesep directoryNames{i} '.csv'], neighbrsData, neighbrHeaders);
    %catch e
        %disp(e.identifier);
        %disp(e.message);
    %end
end

