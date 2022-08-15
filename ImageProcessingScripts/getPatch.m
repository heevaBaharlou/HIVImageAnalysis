function [CellId, patch1, patch2, patch3, patch4] = ...
    getPatch(cellMask1, cellMask2, cellMask3, cellMask4, expansion)
%Load image size
[sr,sc]  = size(cellMask1);

%List of pixel indices (linear)
props     = regionprops(cellMask1,'PixelIdxList');
lenIDs    = unique(cellMask1);
len       = double(lenIDs(lenIDs ~= 0));

%This first part is for percent touching and number neighbor calculation on
%a fixed pixel expansion
CellId = len';
%Get coordinates of each CellId
[r,c] = cellfun(@(x) ind2sub([sr sc],x),{props(CellId).PixelIdxList},'UniformOutput', false);
%The below conditions check the min and max conditions for each.Do NOT CHANGE
rmax = cellfun(@(x) min(sr,max(x) + (expansion)), r);
rmin = cellfun(@(x) max(1,min(x)  - (expansion)),r);
cmax = cellfun(@(x) min(sc,max(x) + (expansion)),c);
cmin = cellfun(@(x) max(1,min(x)  - (expansion)),c);
idxr = cellfun(@(x,y) x:y, num2cell(rmin),num2cell(rmax),'UniformOutput',false);
idxc = cellfun(@(x,y) x:y, num2cell(cmin),num2cell(cmax),'UniformOutput',false);
patch1 = cellfun(@(x,y) cellMask1(x,y),idxr,idxc,'UniformOutput',false);
patch2 = cellfun(@(x,y) cellMask2(x,y),idxr,idxc,'UniformOutput',false);
patch3 = cellfun(@(x,y) cellMask3(x,y),idxr,idxc,'UniformOutput',false);
patch4 = cellfun(@(x,y) cellMask4(x,y),idxr,idxc,'UniformOutput',false);
