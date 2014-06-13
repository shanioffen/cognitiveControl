function saveROItSeries(subj,expt)
% this function has to be run from within a mrLR session
% and you need to have loaded all the ROIs you want


  
gnum = viewGet(getMLRView,'groupNum','Concatenation');
concatInfo = viewGet(getMLRView,'concatInfo',1,gnum);
numRoi = viewGet(getMLRView,'numrois')
roiTimeCourses = tseriesROI(getMLRView,gnum,1:numRoi,1);

switch expt
  case 'detect'
    saveDir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/roiTseries/';
  case 'memory'
    saveDir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/roiTseries/';
end

for iRoi = 1:numRoi

  ROI = viewGet(getMLRView,'roi',iRoi);
  if length(ROI.coords)< 100
    disp(sprintf('ROI %s too small, being eliminated',ROI.name))
    roiTseries = NaN*ones(size(roiTimeCourses{iRoi}));
  else
    roiTseries = roiTimeCourses{iRoi};
  end
  
  save([saveDir subj '_' expt '_' ROI.name],'roiTseries','concatInfo')
  clear roiTseries 

end
