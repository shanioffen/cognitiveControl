function mkTables

tableDir = '/Users/shani/Dropbox/FEF_IPS_paper/tables/raw/';
bootdir = '/Users/shani/NYU/fMRI_data/GLM_output/RoB_ROIanalysis/bootstrap/';
sizeDir = '/Users/shani/NYU/fMRI_data/GLM_output/RoB_ROIanalysis/ROIsize/';



% make tables of bootstrap results
load([bootdir 'bootResults']);

ROIlist = bootResults.ROIlist; % so know what order were run in for boots;
subjList = bootResults.subjList;
exptList = bootResults.exptList;

numRoi = length(ROIlist);
numExpt = length(exptList);
numSubj = length(subjList);

total = length(subjList)*length(exptList);




dlmwrite([tableDir 'FDR_bootPvalCorrect.csv'],bootResults.pGreaterCorr,',');
dlmwrite([tableDir 'FDR_bootPvalWrong.csv'],bootResults.pGreaterWrong,',');
dlmwrite([tableDir 'FDR_dValCorr.csv'],reshape(bootResults.dCorr,total,length(ROIlist)),',');
dlmwrite([tableDir 'FDR_dValWrong.csv'],reshape(bootResults.dWrong,total,length(ROIlist)),',');  
dlmwrite([tableDir 'FDR_bootPvalDif.csv'],squeeze(bootResults.percentGreaterDifCorr),',');  
dlmwrite([tableDir 'FDR_dValDif.csv'],squeeze(bootResults.dDifCorr),',');    

% also get r2 for each subj/expt/ROI and make a table
for iExpt = 1:numExpt
  expt = exptList{iExpt};
  for iSubj = 1:numSubj
    subj = subjList{iSubj};
    for iROI = 1:numRoi
      ROI = ROIlist{iROI};
      glm = loadGLM(subj,expt,ROI);    
      r2(iExpt,iSubj,iROI) = glm.epochR2;
    end
  end
end
r2 = reshape(r2,total,numRoi); % want in same shape and order at boot results;
dlmwrite([tableDir 'R2table.csv'],r2,',');

% also get ROIsizes for each subj/ROI and make a table
keepsize = NaN*ones(numSubj,numRoi);

for iSubj = 2:numSubj
  subj = subjList{iSubj};
  for iRoi = 1:numRoi
    ROI = ROIlist{iRoi};
    
    sizeName = [sizeDir 'ROIsize_' subj '_' ROI '.mat'];
    load(sizeName)
    keepSize(iSubj,iRoi) = roiSize;
    clear roiSize;
    
  end
end
keepSize(1,:) = nanmean(keepSize);
dlmwrite([tableDir 'FDR_ROISize.csv'],keepSize,',')



      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function glm = loadGLM(subj, expt, ROI)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch expt
  case{'detect'}
    glmdir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/roiBetas/';
  case{'memory'}
    glmdir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/roiBetas/';
end

dataName = [glmdir 'glm_' subj '_' expt '_' ROI '_correctTrialS2plus_estHRF.mat'];

load(dataName)
if(exist('glm'))
  glm = glm;
else
  glm = [];
end

