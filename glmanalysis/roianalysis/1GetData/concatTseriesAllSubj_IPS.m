function concatTseriesAllSubj_IPS(expt)

  
subjList = {'JG','DS','LM','RS','SO'};
numSubj = length(subjList);
ROIlist = {...
    'newLeftIPS_restCan25', 'newRightIPS_restCan25',...
    'newLeftIPS_restCan30', 'newRightIPS_restCan30', ...
    'newLeftIPS_restCan35', 'newRightIPS_restCan35' ...
          };
numRoi = length(ROIlist);

switch expt
  case 'detect'
    saveDir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/roiTseries/';
  case 'memory'
    saveDir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/roiTseries/';
end

for iRoi = 1:numRoi
  ROI = ROIlist{iRoi};
  TC = [];

  for iSubj = 1:numSubj
    subj= subjList{iSubj};
    roiDataName = [saveDir subj '_' expt '_' ROI];
    load(roiDataName);
    % get average ROI time course
    roiTC = mean(roiTseries,2);
    
    
    TC = [TC; roiTC];
    cInf{iSubj} = concatInfo;
    clear roiTC concatInfo;
  end % going through subjects
  roiTseries = TC;
  concatInfo = cInf;
  saveDataName = [saveDir 'allSubj_' expt '_' ROI];
  save(saveDataName,'roiTseries','concatInfo');
  clear roiTseries concatInfo
end
