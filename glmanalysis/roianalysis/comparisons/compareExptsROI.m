function compareExptsROI(subj,modelNum,hrfNum)
% compares how well different models match the data in prefrontal ROIs
  
redo = 1; % whether to recalculate the GLMs or just load them

if nargin == 2
  hrfNum = 2;
elseif nargin == 1
  modelNum = 6;
  hrfNum = 2;
end

hrfList = {'Vis','Can','Est'};
numHrf = length(hrfList);
whichHrf = hrfList{hrfNum};

exptList = {'detect','memory'};
numExpt = length(exptList);

roiList = {'leftDMsPCS', 'rightDMsPCS', 'leftDLsPCS', 'rightDLsPCS'};
numRoi = length(roiList);

modelList = {'correctTrial','correctTrialS2plus','s2plusLong','correctTrialS2trip','S2tripLong','none'};  
modelNickname = {'S2single','S2plus','S2plusLong','S2trip','S2tripLong','none'};
numModel = length(modelList);
model = modelList{modelNum};
modnik = modelNickname{modelNum};

figNum = figure;
plotNum = 0;

for iExpt = 1:numExpt
  expt = exptList{iExpt};
  
  switch expt
    case{'detect'}
      datadir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/roiTseries/';
    case{'memory'}
      datadir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/roiTseries/';
  end
  
  % get the start TR, delayDuration, and correctStatus for calculating TTA
  exptTiming = getExptTiming(subj,expt); %subfunction
  
  for iRoi = 1:numRoi
    ROI = roiList{iRoi};
    % load ROI time course and concatInfo
    % this loads two variables: roiTseries and concatInfo
    load([datadir subj '_' expt '_' ROI '.mat']);
    
    % get average ROI time course
    roiTC = mean(roiTseries,2);
    % get relevant concatInfo
    runTransition = concatInfo.runTransition;
    hipassfilter = concatInfo.hipassfilter;
    
    % calculate TTA for plotting:
    dataTTAcorrect = getTTA(roiTC,exptTiming,1); % get TTA for correct
    dataTTAwrong = getTTA(roiTC, exptTiming,-1);
    if ~strcmp(model,'none')
      % must run compareModels to get glm; here just load it
      glm = loadGLM(subj,expt,model,ROI,canHrf);
      
      % now calculate TTA for the model
      modelTTAcorrect = getTTA(glm.model,exptTiming,1);
      modelTTAwrong = getTTA(glm.model,exptTiming,-1);
      
      % now plot
      dpi = glm.betas(2)/(0.5*(glm.betas(1)+glm.betas(3)));
      plotNum = 2*iRoi + (iExpt - 2);
      figure(figNum), subplot(numRoi,numExpt,plotNum);
      errorbar(-4:2:30, dataTTAcorrect(:,1),dataTTAcorrect(:,2),'r.-.');
      hold on;
      legendHandler = errorbar(-4:2:30,modelTTAcorrect(:,1),modelTTAcorrect(:,2),'k.-');
      legendStr = sprintf('dpi %0.2f',dpi);
      %    xlabel('time');
      xlim([-4 30]);
      %    ylabel('MRI signal');
      ylim([-0.2 0.3])
      title(sprintf('corr %s %s %s %0.2f',ROI,expt,modnik,glm.r2));
      legend(legendHandler,legendStr,'Location','NorthWest');
      
    else %plot just data with no model
      % now plot
      plotNum = 2*iRoi + (iExpt - 2);
      figure(figNum), subplot(numRoi,numExpt, plotNum);
      errorbar(-4:2:30, dataTTAcorrect(:,1),dataTTAcorrect(:,2),'r.-.');
      xlim([-4 30]);
      ylim([-0.2 0.3]);
      title(sprintf('%s - %s - %s',subj,expt,ROI));
    
    end
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function glm = loadGLM(subj, expt, model, ROI, whichHrf)
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

load(dataName)
if(exist('glm'))
  glm = glm;
else
  glm = [];
end

  

