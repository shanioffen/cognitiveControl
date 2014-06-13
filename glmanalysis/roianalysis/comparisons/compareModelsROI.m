function compareModelsROI(subj,expt,HRFnum)
% function compareModelsROI(subj,expt,HRFnum)
% enter 1 for Vis HRF, 2 for Canonical HRF, and 3 for estimated HRF
  
% compares how well different models match the data in prefrontal ROIs
  
redo = 1; % whether to recalculate the GLMs or just load them

if nargin==2
  HRFnum = 3;
end

hrfList = {'Vis','Can','Est'};
numHrf = length(hrfList);
whichHRF = hrfList{HRFnum};

roiList = {'leftDMsPCS', 'rightDMsPCS', 'leftDLsPCS', 'rightDLsPCS'};
numRoi = length(roiList);

modelList = {'correctTrial','correctTrialS2plus','s2plusLong','correctTrialS2trip','S2tripLong'};  
modelNickname = {'S2single','S2plus','S2plusLong','S2trip','S2tripLong'};
numModel = length(modelList);

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
  dataTTAcorrect = getTTA(roiTC,exptTiming,1); % get TTA for correct
  dataTTAwrong = getTTA(roiTC, exptTiming,-1);
  for iModel = 2; %1:numModel
    model = modelList{iModel};
    modnik = modelNickname{iModel};
    
    % do GLM on that time course if haven't already
    if(redo)
      % load model DM (for now only look at trials separated by correct/incorrect)
      DMlong = getDM(subj,expt,model); % subfunction
      
      % load the HRF that will be used to convolve the DM
      hrfParams = loadHRFparams(subj,whichHRF,ROI,model); % subfunction;
      hrfParams.tmax = 16;
      hrf = spmHRF_so(2,hrfParams); % separate fxn;
      glm.hrf = hrf; % save it to structure
      glm.hrfParams = hrfParams; %keep the paramters used to make the HRF
      
      %convolve DM by HRF one run at a time and stack
      scm = [];
      % test how long this takes
      for runnum = 1:size(runTransition,1) 
        DM = [];
        % convolve DM with HRF
        DM = DMlong(runTransition(runnum,1):runTransition(runnum,2),:);
        m = convn(DM, hrf(:,1));
        m = m(1:length(DM),:);
        % remove mean 
        m = m-repmat(mean(m), size(m,1), 1);
        % apply the same filter as original data
        m = real(ifft(fft(m) .* repmat(hipassfilter{runnum}', 1, size(m,2)) ));
        %stack
        scm = [scm; m];
      end
      glm.scm = scm; % save the convolved DM
      
      % calculate the GLM and get betas
      % precalculate the normal equation (this dramatically speeds up things)
      precalcmatrix = ((scm'*scm)^-1)*scm';
      % if this don't work then do pinv
      if sum(isnan(precalcmatrix(:))) == length(precalcmatrix(:))
        disp(sprintf('(riGLMandPlot) Using pseudo inverse to invert convolution matrix'));
        precalcmatrix = pinv(scm);
      end

      glm.betas = precalcmatrix*roiTC;
      % calculate error bars, first get sum-of-squares of residual
      % (in percent signal change)
      glm.model = scm * glm.betas;
      sumOfSquaresResidual = sum((roiTC-glm.model).^2);
      % now calculate the sum-of-squares of that error
      % and divide by the degrees of freedom (n-k where n
      % is the number of timepoints in the scan and k is 
      % the number of timepoints in all the estimated hdr)
      S2 = sumOfSquaresResidual/(length(roiTC)-size(scm,2));
      % now distribute that error to each one of the points
      % in the hemodynamic response according to the inverse
      % of the covariance of the stimulus convolution matrix.
      glm.betaError = sqrt(diag(pinv(scm'*scm))*S2);
      % calculate variance accounted for by the estimated hdr
      glm.r2 = (1-sumOfSquaresResidual./sum(roiTC.^2));

      
      % save all that info
      saveGLM(subj,expt,model,ROI,glm,whichHRF) % subfunction ;
    else % if don't redo, then load
      glm = loadGLM(subj,expt,model,ROI,whichHRF);
    end % if(redo)
    
    % now calculate TTA for the model
    modelTTAcorrect = getTTA(glm.model,exptTiming,1);
    modelTTAwrong = getTTA(glm.model,exptTiming,-1);
    
    % get epoch-based R2
    [epochR2 resid] = calcVarAccnt(dataTTAcorrect, modelTTAcorrect);
        
    % now plot
    plotNum = plotNum + 1; %can't remember quite how this goes
    figure(figNum), subplot(numRoi,numModel,plotNum);
    legendHandle(1) = errorbar(-4:2:30, dataTTAcorrect(:,1),dataTTAcorrect(:,2),'r.-.');
    legendStr{1} = sprintf('Data - %s %s',subj,expt);
    hold on;
    legendHandle(2) = errorbar(-4:2:30,modelTTAcorrect(:,1),modelTTAcorrect(:,2),'k.-');
    legendStr{2} = 'Model';
    %    xlabel('time');
    xlim([-4 30]);
    %    ylabel('MRI signal');
    ylim([-0.2 0.3])
    title(sprintf('corr %s %s %0.2f',ROI,modnik,epochR2));
    %    legend(legendHandle,legendStr);
    
    plotMore = 0;
    if(plotMore)
      % plot incorrect
      figure(figNum+1), subplot(numRoi,numModel,plotNum);
      legendHandle(1) = errorbar(-4:2:30, dataTTAwrong(:,1),dataTTAwrong(:,2),'r.-.');
      legendStr{1} = sprintf('Data - %s %s',subj,expt);
      hold on;
      legendHandle(2) = errorbar(-4:2:30,modelTTAwrong(:,1),modelTTAwrong(:,2),'k.-');
      legendStr{2} = 'Model';
      %    xlabel('time');
      xlim([-4 30]);
      %    ylabel('MRI signal');
      ylim([-0.2 0.3])
      title(sprintf('inc %s %s %0.2f',ROI,modnik,glm.r2));
      %    legend(legendHandle,legendStr);
      
      % show corr vs incorr
      figure(figNum+2), subplot(1,numRoi,iRoi);
      legendHandle(1) = errorbar(-4:2:30, dataTTAcorrect(:,1),dataTTAcorrect(:,2),'r-');
      legendStr{1} = sprintf('correct - %s',subj);
      hold on;
      legendHandle(2) = errorbar(-4:2:30, dataTTAwrong(:,1),dataTTAwrong(:,2),'g-');
      legendStr{2} = sprintf('Incorrect -%s',expt);
      %    xlabel('time');
      xlim([-4 30]);
      %    ylabel('MRI signal');
      ylim([-0.2 0.3])
      title(sprintf('corr %s %s %0.2f',ROI,modnik,glm.r2));
      legend(legendHandle,legendStr);
    end
    
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function timing = getExptTiming(subj, expt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

saveDir = '/Users/shani/NYU/fMRI_data/GLM_output/DM/'; %where to save the design matrices
dataName = [saveDir 'exptTiming_' expt '_' subj];
load(dataName);
if(exist('timing'))
  timing = timing;
else
  timing = [];
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
function saveGLM(subj, expt, model, ROI, glm, whichHRF)
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

save(dataName, 'glm');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function glm = loadGLM(subj, expt, model, ROI, whichHRF)
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

