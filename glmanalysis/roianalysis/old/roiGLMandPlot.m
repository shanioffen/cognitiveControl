function roiGLMandPlot(subj,expt)
  
% compares how well different models match the data in prefrontal ROIs
  
redo = 1; % whether to recalculate the GLMs or just load them


roiList = {'leftDMsPCS', 'leftDLsPCS', 'rightDMsPCS', 'rightDLsPCS'};
numRoi = length(roiList);
modelList = {'correctTrial','correctTrialS2plus','correctTrialS2minus'};  
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
  for iModel = 1:numModel
    model = modelList{iModel};
    
    % do GLM on that time course if haven't already
    if(redo)
      % load model DM (for now only look at trials separated by correct/incorrect)
      DMlong = getDM(subj,expt,model); % subfunction
      
      % load the HRF that will be used to convolve the DM
      hrfParams = loadHRFparams(subj); % subfunction;
      hrfParams.tmax = 16;
      hrf = spmHRF_so(2,hrfParams); % separate fxn
      
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
      saveGLM(subj,expt,model,ROI,glm) % subfunction ;
    else % if don't redo, then load
      glm = loadGLM(subj,expt,model,ROI);
    end % if(redo)
    
    % now calculate TTA for the model
    modelTTAcorrect = getTTA(glm.model,exptTiming,1);
    modelTTAwrong = getTTA(glm.model,exptTiming,-1);
    
    % now plot
    plotNum = plotNum + 1; %can't remember quite how this goes
    figure(figNum), subplot(numRoi,numModel,plotNum);
    legendHandle(1) = errorbar(-4:2:30, dataTTAcorrect(:,1),dataTTAcorrect(:,2),'r.-.');
    legendStr{1} = sprintf('Data - %s %s',subj,expt);
    hold on;
    legendHandle(2) = errorbar(-4:2:30,modelTTAcorrect(:,1),modelTTAcorrect(:,2),'k.-');
    legendStr{2} = 'Model';
    xlabel('time');
    xlim([-4 30]);
    ylabel('MRI signal');
    %    ylim([-0.5 1])
    title(sprintf('TTA %s %s',ROI,model));
    legend(legendHandle,legendStr);
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hrfParams = loadHRFparams(subj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hrfDir = '/Users/shani/NYU/fMRI_data/GLM_output/HRF_est/parameters/';

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveGLM(subj, expt, model, ROI, glm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch expt
  case{'detect'}
    glmdir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/roiBetas/';
  case{'memory'}
    glmdir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/roiBetas/';
end

dataName = [glmdir 'glm_' subj '_' expt '_' ROI '_' model];
save(dataName, 'glm');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function glm = loadGLM(subj, expt, model, ROI)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch expt
  case{'detect'}
    glmdir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/roiBetas/';
  case{'memory'}
    glmdir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/roiBetas/';
end

dataName = [glmdir 'glm_' subj '_' expt '_' ROI '_' model];
load(dataName)
if(exist('glm'))
  glm = glm;
else
  glm = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TTA = getTTA(TC,timing,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the TTA for a time course, given the timing of the expt
% can get TTA only for correct trials (right = 1), only incorrect (right = -1), 
% or all (right = 0) Default is 1
  
if nargin < 3
  right = 1;
elseif nargin == 3
  right = varargin{1};
end

% find which trials meet criteria:
correct = timing.correct;
if right == 1
  takeIndx = find(correct == 1);
elseif right == -1
  takeIndx = find(correct == 0);
elseif right == 0
  takeIndx = 1:length(correct);
end
  
% check the delays for those trials
starts = timing.startTR(takeIndx);
delays = timing.delayTR(takeIndx);
cutoff = 6; % only delays 12 sec and longer should be graphed;
keepIndx = find(delays >= cutoff);

starts = starts(keepIndx);

numTrials = length(keepIndx);
epochTseries = NaN*ones(numTrials,18); 
% 18 time points = the longest delay period (8) + ITI  (8) + 2 volumes before trial start


for iTrial = 1:numTrials
  if starts(iTrial) == 1 % for first trial, special case
    startPnt = 1;
    endPnt = startPnt + 15; % since not get 2 pnts before trial start
    epochTseries(iTrial,3:18) = TC(startPnt:endPnt);
  else
    startPnt = starts(iTrial)-2; 
    endPnt = min(startPnt + 17, length(TC)); % don't go off other edge either 
    epochLength = endPnt - startPnt + 1;
    epochTseries(iTrial,1:epochLength) = TC(startPnt:endPnt);
  end
end

% there's a problem with how I'm doing this: when getting the final trial before a scan
% ended, I'm going to be taking data starting in the next scan too when I take all 18 pnts
% might need to go through things scan by scan, but that's a pain. Or figure out a way
% to easily calculate where the transitions happen, by subtracting starts from itself
% shifted one, and then deal with that.
TTA(:,1) = nanmean(epochTseries);
TTA(:,2) = nanstd(epochTseries)./sqrt(numTrials);


  

