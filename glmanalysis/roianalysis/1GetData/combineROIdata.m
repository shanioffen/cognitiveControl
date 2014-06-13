function combineROIdata(restrict,ROIlist,newROIlist)

if nargin == 1
  ROIlist = [];
elseif nargin == 0
  restrict = 1;
  ROIlist = [];
end

  
subjList = {'JG','DS','LM','RS','SO'};
numSubj = length(subjList);
exptList = {'detect','memory'};
numExpt = length(exptList);

if isempty(ROIlist)
ROIlist = {...
           'leftDMsPCS',...
           'rightDMsPCS',...
           'leftDLsPCS',...
           'rightDLsPCS',...
           'leftIPCS',...
           'rightIPCS',...
           'leftAntCS',...
           'rightAntCS',...
           'leftCing',...
           'rightCing',...
           'leftPostIPS'...
           'rightPostIPS'...
          }

newROIlist = {...
           'bilatDMsPCS',...
           'bilatDLsPCS',...
           'bilatIPCS',...
           'bilatAntCS',...
           'bilatCing',...
           'bilatPostIPS'...
          };

end



numRoi = length(ROIlist);

for iExpt = 1:numExpt
  expt = exptList{iExpt};
  switch expt
    case 'detect'
      saveDir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/roiTseries/';
    case 'memory'
      saveDir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/roiTseries/';
  end
  for iRoi = 1:2:numRoi
    ROIlh = ROIlist{iRoi}
    ROIrh = ROIlist{iRoi+1}
    newROI = newROIlist{(iRoi+1)/2}
    allTC = [];
    for iSubj = 1:numSubj
      subj= subjList{iSubj};
      
      roiDataNameLH = [saveDir subj '_' expt '_' ROIlh];
      roiDataNameRH = [saveDir subj '_' expt '_' ROIrh];
      saveDataName = [saveDir subj '_' expt '_' newROI];
      
      if restrict
        roiDataNameLH = [roiDataNameLH '_restCan25'];
        roiDataNameRH = [roiDataNameRH '_restCan25'];        
        saveDataName = [saveDataName '_restCan25'];        
      end
      
      load(roiDataNameLH);
      TC = roiTseries;
      clear roiTseries;
      
      load(roiDataNameRH);
      TC = [TC roiTseries];
      clear roiTseries;
      
      roiTseries = TC;
      save(saveDataName,'roiTseries','concatInfo');

      % also concatenate mean across subjects
      allTC = [allTC; nanmean(TC,2)];
      cInf{iSubj} = concatInfo;
      
      clear roiTseries concatInfo TC
    
    end % going through subjects
    roiTseries = allTC;
    concatInfo = cInf;
    saveAllDataName = [saveDir 'allSubj_' expt '_' newROI];
    if restrict
      saveAllDataName = [saveAllDataName '_restCan25'];
    end
    save(saveAllDataName,'roiTseries','concatInfo');
    clear roiTseries concatInfo
  end % going throughROIs
end % going through expts

  
