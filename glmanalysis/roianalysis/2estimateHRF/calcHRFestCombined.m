function calcHRFestCombined
% need to combine estimates from the 4 ROIs so have one to use for the maps

roiList = {'leftDMsPCS', 'rightDMsPCS', 'leftDLsPCS', 'rightDLsPCS'};
numRoi = length(roiList);

subjList = {'JG','DS','LM','RS','SO'};
numSubj = length(subjList);

% fixing one model but leave all this here just in case...
modelList = {'correctTrial','correctTrialS2plus','s2plusLong','correctTrialS2trip','S2tripLong'};  
modelNickname = {'S2single','S2plus','S2plusLong','S2trip','S2tripLong'};
numModel = length(modelList);
model = modelList{2}; modelnik = modelNickname{2};

glmDir = '/Users/shani/NYU/fMRI_data/GLM_output/RoB_ROIanalysis/';

figNum = figure;
plotNum = 0;

for iSubj = 1:numSubj
  subj = subjList{iSubj};
  for iRoi = 1:numRoi
    ROI = roiList{iRoi};
    
    % load the estimates for each Roi
    hrfParams = loadHRFparams(subj,'Est',ROI,model);
    hrfParams.tmax = 16;
    hrf = spmHRF_so(2,hrfParams); % separate fxn;
    rdelay(iRoi) = hrfParams.rdelay;
    udelay(iRoi) = hrfParams.udelay;
    clear hrfParams;
    
    % plot across ROIs to compare
    plotNum = plotNum + 1;
    figure(figNum), subplot(numSubj, numRoi, plotNum);
    h = plot(hrf,'k.-');
    s = sprintf('resp %.2f under %.2f',rdelay(iRoi),udelay(iRoi)); 
    title(sprintf('subj %s ROI %s',subj,ROI));
    legend(h,s);
    ylim([-0.1 0.4])
    
  end
  % average across ROIs and save as HRF est
  estParams.rdelay = mean(rdelay);
  estParams.udelay = mean(udelay);
  dataName = [glmDir 'HRFest_' subj];
  save(dataName,'estParams');
  clear hrfParams rdelay udelay
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hrfParams = loadHRFparams(subj,whichHRF,ROI,model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hrfDir = '/Users/shani/NYU/fMRI_data/GLM_output/HRF_est/parameters/';

% estimated HRF using detect and memory both, since need to use
% same HRF for both experiments.
glmHRFdir = '/Users/shani/NYU/fMRI_data/GLM_output/RoB_ROIanalysis/';

switch whichHRF
  case 'Can'
    hrfParams.rdelay = 6;
    hrfParams.udelay = 16;
    hrfParams.udispersion = 1;
    hrfParams.rdispersion = 1;
    if strcmp(subj,'allSubj')
      hrfParams.rdelay = 6*ones(5,1);
      hrfParams.udelay = 16*ones(5,1);
      hrfParams.udispersion = 1*ones(5,1);
      hrfParams.rdispersion = 1*ones(5,1);
    end
  
  case 'Vis'
    if strcmp(subj,'allSubj')
      % have to deal separately with each subject
      subjList = {'JG','DS','LM','RS','SO'};%switched order because looking on JG brain so put first
      for iSub = 1:length(subjList)
        sub = subjList{iSub};
        HRFfileName = [hrfDir 'SPM_HRF_params_' sub '_V1V2V3_restrict_21'];
        load(HRFfileName) % this loads the variable bestfit;
        hrfParams.rdelay(iSub,1) = bestfit.params(1);
        hrfParams.udelay(iSub,1) = bestfit.params(2);
        hrfParams.udispersion(iSub,1) = bestfit.params(3);
        hrfParams.rdispersion(iSub,1) = bestfit.params(4);
        clear bestfit;
      end
    else % just load for individual subject if not doing fixed effects
      HRFfileName = [hrfDir 'SPM_HRF_params_' subj '_V1V2V3_restrict_21'];
      load(HRFfileName) % this loads the variable bestfit
      hrfParams.rdelay = bestfit.params(1);
      hrfParams.udelay = bestfit.params(2);
      hrfParams.udispersion= bestfit.params(3);
      hrfParams.rdispersion = bestfit.params(4);
    end
  case 'Est'
    if strcmp(subj,'allSubj')
      % have to deal separately with each subject
      subjList = {'JG','DS','LM','RS','SO'};%switched order because looking on JG brain so put first
      for iSub = 1:length(subjList)
        sub = subjList{iSub};
        dataName = [glmHRFdir 'glm_' subj '_' ROI '_' model '_estHRF_bothExpt'];      
        load(dataName); % this loads variable glm with field glm.hrfParams
        hrfParams.rdelay(iSub,1) = glm.hrfParams.rdelay;
        hrfParams.udelay(iSub,1) = glm.hrfParams.udelay;
        hrfParams.udispersion(iSub,1) = 1;
        hrfParams.rdispersion(iSub,1) = 1;
        clear glm;
      end
    else % just load for individual subject if not doing fixed effects
      dataName = [glmHRFdir 'glm_' subj '_' ROI '_' model '_estHRF_bothExpt'];      
      load(dataName); % this loads variable glm with field glm.hrfParams
      hrfParams = glm.hrfParams;
    end
end
