function concatTseriesAllSubjmoreROI(expt,restrict,ROIlist)

if nargin == 2
  ROIlist = [];
end

  
subjList = {'JG','DS','LM','RS','SO'};
numSubj = length(subjList);

if isempty(ROIlist)
  ROIlist = {'leftAntCS','leftDMsPCS','leftDLsPCS','leftIPCS',...
           'rightAntCS','rightDMsPCS','rightDLsPCS','rightIPCS',...
           'leftCing','leftSupBig','leftInfBig','leftPostIPS'...
           'rightCing','rightSupBig','rightInfBig','rightPostIPS'};
end

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
