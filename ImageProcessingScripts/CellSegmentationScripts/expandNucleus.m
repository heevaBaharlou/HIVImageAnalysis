function cellFinal = expandNucleus(nuclearMask,markerMask)

empty = cell2mat(struct2cell(regionprops(logical(markerMask),nuclearMask,'MaxIntensity')));
pixelIdxList = regionprops(logical(markerMask),'PixelIdxList');

for i = 1:length(pixelIdxList)
    if empty(i) == 0
        markerMask(pixelIdxList(i).PixelIdxList) = 0;
    end
end

combined = nuclearMask+markerMask;
se = strel('disk',10);
combined = imclose(combined,se);

DM1 = bwdist(nuclearMask);
DM1(combined==0) = Inf;

cell1 = watershed(DM1)>0;
cell1(combined==0) = 0;

cellEmpty1 = cell2mat(struct2cell(regionprops(logical(cell1),nuclearMask,'MaxIntensity')));
pixelIdxList = regionprops(logical(cell1),'PixelIdxList');

for i = 1:length(pixelIdxList)
    if cellEmpty1(i) == 0
        cell1(pixelIdxList(i).PixelIdxList) = 0;
    end
end

DM2 = bwdist(cell1);
DM2(combined==0) = Inf;

cell2 = watershed(DM2)>0;
cell2(combined==0) = 0;

cellEmpty2 = cell2mat(struct2cell(regionprops(logical(cell2),nuclearMask,'MaxIntensity')));
pixelIdxList = regionprops(logical(cell2),'PixelIdxList');

for i = 1:length(pixelIdxList)
    if cellEmpty2(i) == 0
        cell2(pixelIdxList(i).PixelIdxList) = 0;
    end
end

cell2 = imfill(cell2,'holes');
cell2 = bwlabel(cell2);

[dist,labels] = bwdist(cell2);
cellBinary = dist<=1;
cellLabel = cell2(labels);
cellFinal = zeros(size(cellLabel));
cellFinal(cellBinary) = cellLabel(cellBinary);
end
