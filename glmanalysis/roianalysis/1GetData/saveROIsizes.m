function saveROIsizes(subj)
% this function has to be run from within a mrLR session
% and you need to have loaded all the ROIs you want

gnum = viewGet(getMLRView,'groupNum','Concatenation');
concatInfo = viewGet(getMLRView,'concatInfo',1,gnum);
numRoi = viewGet(getMLRView,'numrois')

saveDir = '/Users/shani/NYU/fMRI_data/GLM_output/RoB_ROIanalysis/';

for iRoi = 1:numRoi

  ROI = viewGet(getMLRView,'roi',iRoi);
  ROIname = ROI.name
  roiSize = length(ROI.coords)
  save([saveDir 'ROIsize_' subj '_' ROIname],'roiSize')
  clear roiSize

end
