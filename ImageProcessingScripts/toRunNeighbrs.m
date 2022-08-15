workDir = pwd;
files = dir;
directoryNames = {files([files.isdir]).name};
directoryNames = directoryNames(~ismember(directoryNames,{'.','..'}));

toRun = 36;

for i = 1:length(directoryNames)
    disp(i)
    %try
%        [num,txt,raw] = xlsread([workDir filesep directoryNames{i} filesep 'cellData.csv']);
%        [nCell,~] = size(num);

        % Cell masks
        cd3CellMask = double(imread([workDir filesep directoryNames{i} filesep 'cd3CellMask.tif']));
        cd11cCellMask = double(imread([workDir filesep directoryNames{i} filesep 'cd11cCellMask.tif']));
        fxiiiaCellMask = double(imread([workDir filesep directoryNames{i} filesep 'fxiiiaCellMask.tif']));
        otherCellMask = logical(imread([workDir filesep directoryNames{i} filesep 'otherNucMask.tif']));
        otherCellMask = bwlabel(bwpropfilt(otherCellMask,'EquivDiameter',[10 50]));
        
        % weird way to get cell counts
        cd3CellArea = regionprops(cd3CellMask, 'Area');
        cd3CellArea = cell2mat(struct2cell(cd3CellArea));
        cd3Count = length(cd3CellArea);
        cd11cCellArea = regionprops(cd11cCellMask, 'Area');
        cd11cCellArea = cell2mat(struct2cell(cd11cCellArea));
        cd11cCount = length(cd11cCellArea);
        fxiiiaCellArea = regionprops(fxiiiaCellMask, 'Area');
        fxiiiaCellArea = cell2mat(struct2cell(fxiiiaCellArea));
        fxiiiaCount = length(fxiiiaCellArea);
        otherCellArea = regionprops(otherCellMask, 'Area');
        otherCellArea = cell2mat(struct2cell(otherCellArea));
        otherCellCount = length(otherCellArea);
        nCell = cd3Count+cd11cCount+fxiiiaCount+otherCellCount;
        
        clearvars cd3CellArea cd11cCellArea fxiiiaCellArea otherCellArea
        
        %% Neighbours
        expansion = 1; % neighbourhood expansion radius
        se = strel('disk', expansion);
        
        neighbrs = zeros(nCell,100);
        currentMax = 0;
        
        % CD3
        [cd3CellId, cd3Cell_cd3Patch, cd3Cell_cd11cPatch, cd3Cell_fxiiiaPatch, cd3Cell_otherPatch] = ...
            getPatch(cd3CellMask, cd11cCellMask, fxiiiaCellMask, otherCellMask, expansion);
        cd3Extended = cellfun(@(x,y) imdilate(x==y,se,'same'),cd3Cell_cd3Patch,num2cell(cd3CellId),'UniformOutput',false);
        
        cd3_cd3Overlap = cellfun(@(x,y) x(y),cd3Cell_cd3Patch,cd3Extended,'UniformOutput',false);
        cd3_cd3Neighbrs = cellfun(@(x) unique(x),cd3_cd3Overlap,'UniformOutput',false);
        cd3_cd11cOverlap = cellfun(@(x,y) x(y),cd3Cell_cd11cPatch,cd3Extended,'UniformOutput',false);
        cd3_cd11cNeighbrs = cellfun(@(x) unique(x),cd3_cd11cOverlap,'UniformOutput',false);
        cd3_fxiiiaOverlap = cellfun(@(x,y) x(y),cd3Cell_fxiiiaPatch,cd3Extended,'UniformOutput',false);
        cd3_fxiiiaNeighbrs = cellfun(@(x) unique(x),cd3_fxiiiaOverlap,'UniformOutput',false);
        cd3_otherOverlap = cellfun(@(x,y) x(y),cd3Cell_otherPatch,cd3Extended,'UniformOutput',false);
        cd3_otherNeighbrs = cellfun(@(x) unique(x),cd3_otherOverlap,'UniformOutput',false);
        
        for j = 1:cd3Count
            neighbrs1 = cd3_cd3Neighbrs{j}';
            neighbrs1 = neighbrs1(neighbrs1~=0 & neighbrs1~=j);
            
            neighbrs2 = cd3_cd11cNeighbrs{j}';
            neighbrs2 = neighbrs2(neighbrs2~=0 & neighbrs2~=j) + cd3Count;
            
            neighbrs3 = cd3_fxiiiaNeighbrs{j}';
            neighbrs3 = neighbrs3(neighbrs3~=0 & neighbrs3~=j) + cd3Count + cd11cCount;
            
            neighbrs4 = cd3_otherNeighbrs{j}';
            neighbrs4 = neighbrs4(neighbrs4~=0 & neighbrs4~=j) + cd3Count + cd11cCount + fxiiiaCount;
            
            allNeighbrs = sort([neighbrs1 neighbrs2 neighbrs3 neighbrs4]);
            
            if length(allNeighbrs) > currentMax
                currentMax = length(allNeighbrs);
            end
            
            neighbrs(j,1:length(allNeighbrs)) = allNeighbrs;
        end
        
        %% Free up memory
        clearvars cd3_cd3Overlap cd3_cd11cOverlap cd3_fxiiiaOverlap cd3_otherOverlap
        clearvars cd3_cd3Neighbrs cd3_cd11cNeighbrs cd3_fxiiiaNeighbrs cd3_otherNeighbrs
        clearvars cd3Cell_cd3Patch cd3Cell_cd11cPatch cd3Cell_fxiiiaPatch cd3Cell_otherPatch cd3Extended
        
        % CD11c
        [cd11cCellId, cd11cCell_cd11cPatch, cd11cCell_fxiiiaPatch, cd11cCell_otherPatch, cd11cCell_cd3Patch] = ...
            getPatch(cd11cCellMask, fxiiiaCellMask, otherCellMask, cd3CellMask, expansion);
        cd11cExtended = cellfun(@(x,y) imdilate(x==y,se,'same'),cd11cCell_cd11cPatch,num2cell(cd11cCellId),'UniformOutput',false);
        
        cd11c_cd3Overlap = cellfun(@(x,y) x(y),cd11cCell_cd3Patch,cd11cExtended,'UniformOutput',false);
        cd11c_cd3Neighbrs = cellfun(@(x) unique(x),cd11c_cd3Overlap,'UniformOutput',false);
        cd11c_cd11cOverlap = cellfun(@(x,y) x(y),cd11cCell_cd11cPatch,cd11cExtended,'UniformOutput',false);
        cd11c_cd11cNeighbrs = cellfun(@(x) unique(x),cd11c_cd11cOverlap,'UniformOutput',false);
        cd11c_fxiiiaOverlap = cellfun(@(x,y) x(y),cd11cCell_fxiiiaPatch,cd11cExtended,'UniformOutput',false);
        cd11c_fxiiiaNeighbrs = cellfun(@(x) unique(x),cd11c_fxiiiaOverlap,'UniformOutput',false);
        cd11c_otherOverlap = cellfun(@(x,y) x(y),cd11cCell_otherPatch,cd11cExtended,'UniformOutput',false);
        cd11c_otherNeighbrs = cellfun(@(x) unique(x),cd11c_otherOverlap,'UniformOutput',false);
        
        for j = 1:cd11cCount
            neighbrs1 = cd11c_cd3Neighbrs{j}';
            neighbrs1 = neighbrs1(neighbrs1~=0 & neighbrs1~=j);
            
            neighbrs2 = cd11c_cd11cNeighbrs{j}';
            neighbrs2 = neighbrs2(neighbrs2~=0 & neighbrs2~=j) + cd3Count;
            
            neighbrs3 = cd11c_fxiiiaNeighbrs{j}';
            neighbrs3 = neighbrs3(neighbrs3~=0 & neighbrs3~=j) + cd3Count + cd11cCount;
            
            neighbrs4 = cd11c_otherNeighbrs{j}';
            neighbrs4 = neighbrs4(neighbrs4~=0 & neighbrs4~=j) + cd3Count + cd11cCount + fxiiiaCount;
            
            allNeighbrs = sort([neighbrs1 neighbrs2 neighbrs3 neighbrs4]);
            
            if length(allNeighbrs) > currentMax
                currentMax = length(allNeighbrs);
            end
            
            neighbrs(j+cd3Count,1:length(allNeighbrs)) = allNeighbrs;
        end
        
        %% Free up memory
        clearvars cd11c_cd3Overlap cd11c_cd11cOverlap cd11c_fxiiiaOverlap cd11c_otherOverlap
        clearvars cd11c_cd3Neighbrs cd11c_cd11cNeighbrs cd11c_fxiiiaNeighbrs cd11c_otherNeighbrs
        clearvars cd11cCell_cd3Patch cd11cCell_cd11cPatch cd11cCell_fxiiiaPatch cd11cCell_otherPatch cd11cExtended
        
        % FXIIIa
        [fxiiiaCellId, fxiiiaCell_fxiiiaPatch, fxiiiaCell_otherPatch, fxiiiaCell_cd3Patch, fxiiiaCell_cd11cPatch] = ...
            getPatch(fxiiiaCellMask, otherCellMask, cd3CellMask, cd11cCellMask, expansion);
        fxiiiaExtended = cellfun(@(x,y) imdilate(x==y,se,'same'),fxiiiaCell_fxiiiaPatch,num2cell(fxiiiaCellId),'UniformOutput',false);
        
        fxiiia_cd3Overlap = cellfun(@(x,y) x(y),fxiiiaCell_cd3Patch,fxiiiaExtended,'UniformOutput',false);
        fxiiia_cd3Neighbrs = cellfun(@(x) unique(x),fxiiia_cd3Overlap,'UniformOutput',false);
        fxiiia_cd11cOverlap = cellfun(@(x,y) x(y),fxiiiaCell_cd11cPatch,fxiiiaExtended,'UniformOutput',false);
        fxiiia_cd11cNeighbrs = cellfun(@(x) unique(x),fxiiia_cd11cOverlap,'UniformOutput',false);
        fxiiia_fxiiiaOverlap = cellfun(@(x,y) x(y),fxiiiaCell_fxiiiaPatch,fxiiiaExtended,'UniformOutput',false);
        fxiiia_fxiiiaNeighbrs = cellfun(@(x) unique(x),fxiiia_fxiiiaOverlap,'UniformOutput',false);
        fxiiia_otherOverlap = cellfun(@(x,y) x(y),fxiiiaCell_otherPatch,fxiiiaExtended,'UniformOutput',false);
        fxiiia_otherNeighbrs = cellfun(@(x) unique(x),fxiiia_otherOverlap,'UniformOutput',false);
        
        for j = 1:fxiiiaCount
            neighbrs1 = fxiiia_cd3Neighbrs{j}';
            neighbrs1 = neighbrs1(neighbrs1~=0 & neighbrs1~=j);
            
            neighbrs2 = fxiiia_cd11cNeighbrs{j}';
            neighbrs2 = neighbrs2(neighbrs2~=0 & neighbrs2~=j) + cd3Count;
            
            neighbrs3 = fxiiia_fxiiiaNeighbrs{j}';
            neighbrs3 = neighbrs3(neighbrs3~=0 & neighbrs3~=j) + cd3Count + cd11cCount;
            
            neighbrs4 = fxiiia_otherNeighbrs{j}';
            neighbrs4 = neighbrs4(neighbrs4~=0 & neighbrs4~=j) + cd3Count + cd11cCount + fxiiiaCount;
            
            allNeighbrs = sort([neighbrs1 neighbrs2 neighbrs3 neighbrs4]);
            
            if length(allNeighbrs) > currentMax
                currentMax = length(allNeighbrs);
            end
            
            neighbrs(j+cd3Count+cd11cCount,1:length(allNeighbrs)) = allNeighbrs;
        end
        
        %% Free up memory
        clearvars fxiiia_cd3Overlap fxiiia_cd11cOverlap fxiiia_fxiiiaOverlap fxiiia_otherOverlap
        clearvars fxiiia_cd3Neighbrs fxiiia_cd11cNeighbrs fxiiia_fxiiiaNeighbrs fxiiia_otherNeighbrs
        clearvars fxiiiaCell_cd11cPatch fxiiiaCell_fxiiiaPatch fxiiiaCell_cd3Patch fxiiiaCell_otherPatch fxiiiaExtended

        % Other
        [otherCellId, otherCell_otherPatch, otherCell_cd3Patch, otherCell_cd11cPatch, otherCell_fxiiiaPatch] = ...
            getPatch(otherCellMask, cd3CellMask, cd11cCellMask, fxiiiaCellMask, expansion);
        otherExtended = cellfun(@(x,y) imdilate(x==y,se,'same'),otherCell_otherPatch,num2cell(otherCellId),'UniformOutput',false);
        
        other_cd3Overlap = cellfun(@(x,y) x(y),otherCell_cd3Patch,otherExtended,'UniformOutput',false);
        other_cd3Neighbrs = cellfun(@(x) unique(x),other_cd3Overlap,'UniformOutput',false);
        other_cd11cOverlap = cellfun(@(x,y) x(y),otherCell_cd11cPatch,otherExtended,'UniformOutput',false);
        other_cd11cNeighbrs = cellfun(@(x) unique(x),other_cd11cOverlap,'UniformOutput',false);
        other_fxiiiaOverlap = cellfun(@(x,y) x(y),otherCell_fxiiiaPatch,otherExtended,'UniformOutput',false);
        other_fxiiiaNeighbrs = cellfun(@(x) unique(x),other_fxiiiaOverlap,'UniformOutput',false);
        other_otherOverlap = cellfun(@(x,y) x(y),otherCell_otherPatch,otherExtended,'UniformOutput',false);
        other_otherNeighbrs = cellfun(@(x) unique(x),other_otherOverlap,'UniformOutput',false);
        
        for j = 1:otherCellCount
            neighbrs1 = other_cd3Neighbrs{j}';
            neighbrs1 = neighbrs1(neighbrs1~=0 & neighbrs1~=j);
            
            neighbrs2 = other_cd11cNeighbrs{j}';
            neighbrs2 = neighbrs2(neighbrs2~=0 & neighbrs2~=j) + cd3Count;
            
            neighbrs3 = other_fxiiiaNeighbrs{j}';
            neighbrs3 = neighbrs3(neighbrs3~=0 & neighbrs3~=j) + cd3Count + cd11cCount;
            
            neighbrs4 = other_otherNeighbrs{j}';
            neighbrs4 = neighbrs4(neighbrs4~=0 & neighbrs4~=j) + cd3Count + cd11cCount + fxiiiaCount;
            
            allNeighbrs = sort([neighbrs1 neighbrs2 neighbrs3 neighbrs4]);
            
            if length(allNeighbrs) > currentMax
                currentMax = length(allNeighbrs);
            end
            
            neighbrs(j+cd3Count+cd11cCount+fxiiiaCount,1:length(allNeighbrs)) = allNeighbrs;
        end
        
        %% Free up memory
        clearvars other_cd3Overlap other_cd11cOverlap other_fxiiiaOverlap other_otherOverlap
        clearvars other_cd3Neighbrs other_cd11cNeighbrs other_fxiiiaNeighbrs other_otherNeighbrs
        clearvars otherCell_cd3Patch otherCell_cd11cPatch otherCell_fxiiiaPatch otherCell_otherPatch otherExtended
         
        %% csvGeneration
        cellIdx = 1:nCell;
        
        neighbrHeaders =  cellstr(num2str((0:currentMax)'))';
        neighbrHeaders = strrep(neighbrHeaders,' ','');
        neighbrHeaders = strcat('N', neighbrHeaders);
        neighbrHeaders{1} = 'CellId';
        neighbrsData = [cellIdx' neighbrs(:,1:currentMax)];
        
        csvwrite_with_headers([workDir filesep directoryNames{i} filesep 'neighbrs.csv'], neighbrsData, neighbrHeaders);
        csvwrite_with_headers([workDir filesep 'neighbrs' filesep directoryNames{i} '.csv'], neighbrsData, neighbrHeaders);
    %catch e
        %disp(e.identifier);
        %disp(e.message);
    %end
end

