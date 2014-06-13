function compareHRFs(subj,expt,modelNum)
% compare the effect of using (1) visually measured HRF, (2) canonical HRF,
% and (3) HRF estimated from ROI data along with the betas.
% Because (3) requires a different analysis (lsqnonlin), need to first
% run compareModels.m and compareMod_estHRF.m to calculate everything:
% then this code will just load the glm data and make plots
  
if nargin ==2
  modelNum =2;
end

roiList = {'leftAntCS','leftDLsPCS','leftDMsPCS','leftIPCS','leftPostIPS',...
           'rightAntCS','rightDLsPCS','rightDMsPCS','rightIPCS','rightPostIPS'};
numRoi = length(roiList);
modelList = {'correctTrial','correctTrialS2plus','s2plusLong','correctTrialS2trip','S2tripLong'};  
modelNickname = {'S2-1','S2-2','S2-2L','S2-3','S2-3L'};
numModel = length(modelList);
whichModel = modelList{modelNum};
modnik = modelNickname{modelNum};
hrfList = {'Vis','Can','Est'};
numHrf = length(hrfList);

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

  % calculate TTA for plotting:
  dataTTAcorrect = getTTA(roiTC,exptTiming,1); % get TTA for correct
  dataTTAwrong = getTTA(roiTC, exptTiming,-1);

  for iHrf = 1:numHrf
    whichHRF = hrfList{iHrf};
    % load glm data
    glm = loadGLM(subj,expt,whichModel,ROI,whichHRF);
    betas = glm.betas;
    hrf = glm.hrf;
    r2 = glm.r2;
    model = glm.model;
    
    % keep track of hrfs to compare later
    if iHrf == 1
      hrfVis = hrf;
    elseif iHrf == 2;
      hrfCan = hrf;
    elseif iHrf == 3
      hrfEst = hrf;
    end
    
    % now calculate TTA for the model
    modelTTAcorrect = getTTA(model,exptTiming,1);
    modelTTAwrong = getTTA(model,exptTiming,-1);

    % now plot
    figure(figNum), subplot(numHrf, numRoi, iRoi + (numRoi)*(iHrf-1));
    errorbar(-4:2:30, dataTTAcorrect(:,1),dataTTAcorrect(:,2),'r.-.');
    hold on;
    legendHandle(1) = errorbar(-4:2:30,modelTTAcorrect(:,1),modelTTAcorrect(:,2),'k.-');
    legendStr{1} = sprintf('%.2f | %.2f | %.2f',betas(1),betas(2),betas(3));
    xlim([-4 30]);
    ylim([-0.2 0.45]);
  set(gca,'xticklabel','off','yticklabel','off','xtick',[],'ytick',[])  

    %    title(sprintf('r2=%.2f %s %s %s',r2,ROI,modnik,whichHRF));
    %    legend(legendHandle,legendStr);
    clear legendHandle legendStr
  end % go through HRFs
  
  plotHRFs = 0;
  if(plotHRFs)
    figure(figNum), subplot(numHrf+1,numRoi,iRoi + numHrf*numRoi);
    legendHandle(1) = plot(hrfVis,'r');
    hold on
    legendHandle(2) = plot(hrfCan,'b');
    legendHandle(3) = plot(hrfEst,'k');
    legendStr{1} = 'Vis';
    legendStr{2} = 'Can';
    legendStr{3} = 'Est';
    %  legend(legendHandle,legendStr);  
    %  title(sprintf('HRF est %s - %s',subj,expt));
    set(gca,'xticklabel','off','yticklabel','off','xtick',[],'ytick',[])  
    clear legendHandle legendStr;
  end
end% go through rois



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



  

