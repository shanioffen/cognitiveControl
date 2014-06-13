function compareDurations(subj,expt, modelNum, canHrf)
% compares how well different models match the data in prefrontal ROIs
% NB: MUST RUN compareMode_estHRF FIRST TO CALCULATE; THIS CODE JUST LOADS  


if nargin==3
  canHrf = 0;
end

roiList = {'leftDMsPCS', 'rightDMsPCS', 'leftDLsPCS', 'rightDLsPCS'};
numRoi = length(roiList);
modelList = {'correctTrial','correctTrialS2plus','s2plusLong','correctTrialS2trip','S2tripLong'};  
modelNickname = {'S2single','S2plus','S2plusLong','S2trip','S2tripLong'};
numModel = length(modelList);
numBins = 4;

switch expt
  case{'detect'}
    datadir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/roiTseries/';
  case{'memory'}
    datadir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/roiTseries/';
end

% get the start TR, delayDuration, and correctStatus for calculating TTA
exptTiming = getExptTiming(subj,expt); %subfunction

figNum = figure;
plotNum = 0;
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
  dataTTAcorrect = getTTA(roiTC,exptTiming,1,numBins); % get TTA for correct
  dataTTAwrong = getTTA(roiTC, exptTiming,-1,numBins);
  for iModel = modelNum
    model = modelList{iModel};
    modnik = modelNickname{iModel};
    
    % must run compareModels to calculate glm; this code just loads is
    glm = loadGLM(subj,expt,model,ROI,canHrf);
    
    % now calculate TTA for the model
    modelTTAcorrect = getTTA(glm.model,exptTiming,1,numBins);
    modelTTAwrong = getTTA(glm.model,exptTiming,-1,numBins);
    
    % now plot
    for iBin = 1:numBins
      plotNum = plotNum + 1; %can't remember quite how this goes
      figure(figNum), subplot(numRoi,numBins+1,plotNum);
      vals = 2*(iBin-1)+1;
      errs = 2*iBin;
      legendHandle(1) = errorbar(-4:2:30, dataTTAcorrect(:,vals),dataTTAcorrect(:,errs),'r.-.');
      legendStr{1} = sprintf('Data - %s %s',subj,expt);
      hold on;
      legendHandle(2) = errorbar(-4:2:30,modelTTAcorrect(:,vals),modelTTAcorrect(:,errs),'k.-');
      legendStr{2} = 'Model';
      %    xlabel('time');
      xlim([-4 30]);
      %    ylabel('MRI signal');
      ylim([-0.2 0.4])
      title(sprintf('corr %s %s %0.2f',ROI,modnik,glm.r2));
      %    legend(legendHandle,legendStr);
    end % going through bins
    plotNum = plotNum + 1; % also plot all together;
    figure(figNum), subplot(numRoi,numBins+1,plotNum);
    errorPlot(-4:2:30, dataTTAcorrect(:,1), dataTTAcorrect(:,2),[1 0 0], 0.2*[1 1 1]);
    hold on;
    errorPlot(-4:2:30, dataTTAcorrect(:,3), dataTTAcorrect(:,4),[0 1 0], 0.2*[1 1 1]);
    errorPlot(-4:2:30, dataTTAcorrect(:,5), dataTTAcorrect(:,6),[0 0 1], 0.2*[1 1 1]);
    errorPlot(-4:2:30, dataTTAcorrect(:,7), dataTTAcorrect(:,8),[1 1 0], 1*[1 1 1]);  
    xlim([-4 30])
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
function glm = loadGLM(subj, expt, model, ROI, canHrf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch expt
  case{'detect'}
    glmdir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/roiBetas/';
  case{'memory'}
    glmdir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/roiBetas/';
end

dataName = [glmdir 'glm_' subj '_' expt '_' ROI '_' model];
if canHrf
  dataName = [dataName '_cannonHRF'];
end
load(dataName)
if(exist('glm'))
  glm = glm;
else
  glm = [];
end


  

