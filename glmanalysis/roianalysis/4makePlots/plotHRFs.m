function [epochR2 betasROI] = plotHRFs
% function estimateBetas(subj,expt,[restrict],[HRFnum],[redo],[doplot])
% calculate betas for given Subj,Expt for all ROIs
% plots time courses
% Model chosen is S2plus (modeling the end of hte trial with 2TRs)
%
% can use canonical HRF or estimated HRF.
% for canonical, enter HRFnum==1 or leave blank;
% for estimated, enter HRFnum==2.
% to use estimated HRF, need to first estimate it using code 'estimateHRFs.m'
  
restrict = 1;
doplot =1;

hrfList = {'Can','Est','Vis'};
subjList = {'JG','DS','LM','RS','SO'};
numSubj = length(subjList);

graphDir = '/Users/shani/NYU/NYUdocuments/thesis/figures/Ch4/data/';

figNum = figure;

plotNum = 0;

for iSubj = 1:numSubj
  subj = subjList{iSubj};
  
  for iHrf = 1:length(hrfList)
    whichHRF = hrfList{iHrf};
    
    % load the HRF that will be used to convolve the DM
    hrfParams = loadHRFparams(subj,whichHRF,restrict);
    hrfParams.tmax = 30;
    hrf(:,iHrf) = spmHRF_so(2,hrfParams); % separate fxn;
  end
  
  plotNum = plotNum + 1; % goes across, then down;
  
  figure(figNum), subplot(3,2,plotNum);
  legendHandle(1) = plot(0:2:30,hrf(:,2),'k');
  legendStr{1} = sprintf('Estimated');
  hold on;
  legendHandle(2) = plot(0:2:30,hrf(:,1),'r');
  legendStr{2} = 'Canonical';
  legendHandle(3) = plot(0:2:30,hrf(:,3),'b');
  legendStr{3} = 'Visual';
  title(sprintf('%s',subj));
  ylim([-.2 .8]);
  xlim([0 30]);

  if plotNum ==1
    legend(legendHandle,legendStr);
  end
  
  if iSubj == numSubj
    graphName = [graphDir 'HRFall_' subj];
    print(gcf,'-painters','-dill',graphName);
  end
  
  figure(figNum+1), subplot(3,2,plotNum);
  plot(0:2:30,hrf(:,2),'k');
  title(sprintf('%s',subj))
  
  if iSubj == numSubj
    graphName = [graphDir 'HRFest_' subj];
    print(gcf,'-painters','-dill',graphName);
  end
  
    
end %going through subj



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [runTransition hipassfilter] = getConcatInfo(subj,expt,concatInfo);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch subj
  case 'allSubj'
    runTransition = [];
    hipassfilter = concatInfo{1}.hipassfilter; % same for all;
    for iSubj = 1:length(concatInfo)
      numRuns(iSubj) = length(concatInfo{iSubj}.runTransition);
      runTemp = concatInfo{iSubj}.runTransition;
      if iSubj>1
        updateVal = sum(numRuns(1:iSubj-1)) * 120;
      else
        updateVal = 0;
      end
      runTransition = [runTransition; runTemp+updateVal];
    end
  otherwise
    runTransition = concatInfo.runTransition;
    hipassfilter = concatInfo.hipassfilter;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hrfParams = loadHRFparams(subj,whichHRF,restrict)
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
        dataName = [glmHRFdir 'HRFest_' sub];    
        if restrict
          dataName = [dataName '_restCan25'];
        end
        load(dataName); % this loads variable estParams
        hrfParams.rdelay(iSub,1) = estParams.rdelay;
        hrfParams.udelay(iSub,1) = estParams.udelay;
        hrfParams.udispersion(iSub,1) = 1;
        hrfParams.rdispersion(iSub,1) = 1;
        clear glm;
      end
    else % just load for individual subject if not doing fixed effects
      dataName = [glmHRFdir 'HRFest_' subj];
      if restrict
        dataName = [dataName '_restCan25'];
      end
      load(dataName); % this loads variable estParams 
      hrfParams.rdelay = estParams.rdelay;
      hrfParams.udelay = estParams.udelay;
      hrfParams.udispersion = 1;
      hrfParams.rdispersion = 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DM = getDM(subj, expt, model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

saveDir = '/Users/shani/NYU/fMRI_data/GLM_output/DM/'; %where to save the design matrices
dataName = [saveDir 'DM_' model '_' expt '_' subj];
load(dataName);
if(exist('DM'))
  DM = DM;
else
  DM = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveGLM(subj, expt, model, ROI, glm, whichHRF,restrict)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch expt
  case{'detect'}
    glmdir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/roiBetas/';
  case{'memory'}
    glmdir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/roiBetas/';
end

dataName = [glmdir 'glm_' subj '_' expt '_' ROI '_' model];

switch whichHRF
  case 'Can'
    dataName = [dataName '_cannonHRF'];
  case 'Est'
    dataName = [dataName '_estHRF'];    
end

if restrict
  dataName = [dataName '_rest'];
end


save(dataName, 'glm');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function glm = loadGLM(subj, expt, model, ROI, whichHRF,restrict)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch expt
  case{'detect'}
    glmdir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/roiBetas/';
  case{'memory'}
    glmdir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/roiBetas/';
end

dataName = [glmdir 'glm_' subj '_' expt '_' ROI '_' model];
switch whichHRF
  case 'Can'
    dataName = [dataName '_cannonHRF'];
  case 'Est'
    dataName = [dataName '_estHRF'];    
end

if restrict
  dataName = [dataName '_rest'];
end

load(dataName)
if(exist('glm'))
  glm = glm;
else
  glm = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function epochR2 = getEpochR2(data, model, exptTiming);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modelTTAmtx = getTTA(model,exptTiming,1,4);
dataTTAmtx = getTTA(data,exptTiming,1,4);

modelTTAall = [];
dataTTAall = [];
for iBin = 1:4
  modelTTAall = [modelTTAall; modelTTAmtx(:,2*(iBin-1)+1)];
  dataTTAall = [dataTTAall; dataTTAmtx(:,2*(iBin-1)+1)];
end

% get rid of the NaNs
notNanInd = find(~isnan(dataTTAall));
modelTTA = modelTTAall(notNanInd);
dataTTA = dataTTAall(notNanInd);

[epochR2 resid] = calcVarAccnt(dataTTA, modelTTA);
