function compareSubjects(expt,modelNum,hrfNum)
% for each expt, compare the data across subjects
  
if nargin == 2
  hrfNum = 2;
elseif nargin == 1
  modelNum = 6;
  hrfNum = 2;
end

subjList = {'JG','DS','LM','RS','SO'};
numSubj = length(subjList);

roiList = {'leftDMsPCS', 'rightDMsPCS', 'leftDLsPCS', 'rightDLsPCS'};
numRoi = length(roiList);

modelList = {'correctTrial','correctTrialS2plus','s2plusLong','correctTrialS2trip','S2tripLong','none'};  
modelNickname = {'S2-1','S2-2','S2-2L','S2-3','S2-3L','noMod'};
numModel = length(modelList);
whichModel = modelList{modelNum};
modnik = modelNickname{modelNum};

hrfList = {'Vis','Can','Est'};
numHrf = length(hrfList);
whichHrf = hrfList{hrfNum};

datadir = getdatadir(expt);

figNum = figure;
plotNum = 0;
for iRoi = 1:numRoi
  ROI = roiList{iRoi};
  for iSubj = 1:numSubj
    subj= subjList{iSubj};

    % load ROI time course and concatInfo
    % this loads two variables: roiTseries and concatInfo
    load([datadir subj '_' expt '_' ROI '.mat']);
    % get average ROI time course
    roiTC = mean(roiTseries,2);
    
    % get the start TR, delayDuration, and correctStatus for calculating TTA
    exptTiming = getExptTiming(subj,expt); %subfunction
    % calculate TTA for plotting:
    dataTTAcorrect = getTTA(roiTC,exptTiming,1); % get TTA for correct
    dataTTAwrong = getTTA(roiTC, exptTiming,-1);
    
    % load glm data if plotting model
    if ~strcmp(whichModel,'none')
      glm = loadGLM(subj,expt,whichModel,ROI,whichHrf);
      betas = glm.betas;
      hrf = glm.hrf;
      r2 = glm.r2;
      model = glm.model;
    
      % now calculate TTA for the model
      modelTTAcorrect = getTTA(model,exptTiming,1);
      modelTTAwrong = getTTA(model,exptTiming,-1);
    
      % now plot
      plotNum = plotNum + 1; % goes across and then down
      figure(figNum), subplot(numRoi,numSubj, plotNum);
      errorbar(-4:2:30, dataTTAcorrect(:,1),dataTTAcorrect(:,2),'r.-.');
      hold on;
      legendHandle(1) = errorbar(-4:2:30,modelTTAcorrect(:,1),modelTTAcorrect(:,2),'k.-');
      legendStr{1} = sprintf('%.2f-%.2f-%.2f',betas(1),betas(2),betas(3));
      xlim([-4 30]);
      ylim([-0.2 0.3]);
      title(sprintf('%s r2=%.2f %s %s %s',subj,r2,ROI,modnik,whichHrf));
      legend(legendHandle,legendStr);
      clear legendHandle legendStr
    else %plot just data with no model
      % now plot
      plotNum = plotNum + 1; % goes across and then down
      figure(figNum), subplot(numRoi,numSubj, plotNum);
      errorbar(-4:2:30, dataTTAcorrect(:,1),dataTTAcorrect(:,2),'r.-.');
      xlim([-4 30]);
      ylim([-0.2 0.3]);
      title(sprintf('%s - %s - %s',subj,expt,ROI));
    end %check if model = none
  end % going through subjects
end % going through ROIs




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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function glm = loadGLM(subj, expt, model, ROI, whichHRF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

load(dataName)
if(exist('glm'))
  glm = glm;
else
  glm = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function datadir = getdatadir(expt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch expt
  case{'detect'}
    datadir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/roiTseries/';
  case{'memory'}
    datadir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/roiTseries/';
end



  

