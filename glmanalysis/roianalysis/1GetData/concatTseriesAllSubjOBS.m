function concatTseriesAllSubj(expt,restrict)
  
subjList = {'JG','DS','LM','RS','SO'};
numSubj = length(subjList);

ROItags = {'leftAntCS','leftDLsPCS','leftDMsPCS','leftIPCS','leftPostIPS',...
           'rightAntCS','rightDLsPCS','rightDMsPCS','rightIPCS','rightPostIPS'};
numRoi = length(ROItags);

switch expt
  case 'detect'
    saveDir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/roiTseries/';
  case 'memory'
    saveDir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/roiTseries/';
end

for iRoi = 1:numRoi
  ROI = ROItags{iRoi};
  TC = [];

  for iSubj = 1:numSubj
    subj= subjList{iSubj};
    roiDataName = [saveDir subj '_' expt '_' ROI];
    if restrict
      roiDataName = [roiDataName '_restCan25'];
    end
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
  if restrict
    saveDataName = [saveDataName '_restCan25'];
  end
  
  save(saveDataName,'roiTseries','concatInfo');
  clear roiTseries concatInfo
end
