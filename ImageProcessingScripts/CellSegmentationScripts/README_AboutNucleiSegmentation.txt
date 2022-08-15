For this new batch, due to deconvolution of DAPI the nuclei segmentation was performed differently to the non-deconvolved data acquired previously. 

In matlab or Image we generated the nuclearMask.tif as follows:
- Gaussian blur with 1.5pixel sigma
- Threshold using manual threshold value from thresholds.csv
- Create mask
- Watershed
- Fill holes
- Output 8bit mask called nuclearMask.tif

Old non-deconvolved data was processed using nucleiSegment.m script.