function compareMod_estHRF(subj,expt)
%OBSOLETE: now estimate HRF for detect and memory simultaneously
% since can't use dif HRF for the two expts
% This function will use the data to estimate the HRF instead of using
% the HRF measured in visual cortex.
  
redo = 1; % whether to recalculate the GLMs or just load them

roiList = {'leftDMsPCS', 'rightDMsPCS', 'leftDLsPCS', 'rightDLsPCS'};
numRoi = length(roiList);
modelList = {'correctTrial','correctTrialS2plus','s2plusLong','correctTrialS2trip','S2tripLong'};  
modelNickname = {'S2single','S2plus','S2plusLong','S2trip','S2tripLong'};
numModel = length(modelList);

TR = 2;
numBins = 1; % not plotting all durations

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
  dataTTAcorrect = getTTA(roiTC,exptTiming,1,numBins); % get TTA for correct
  dataTTAwrong = getTTA(roiTC, exptTiming,-1,numBins);
  for iModel = 1:numModel
    model = modelList{iModel};
    modnik = modelNickname{iModel};
    
    % do GLM on that time course if haven't already
    if(redo)
      % load model DM (for now only look at trials separated by correct/incorrect)
      DM = getDM(subj,expt,model); % subfunction
      
      % initialize variables
      [initialVals lb ub] = setInitialVals(DM); % subfunction
      
      % do lsqnonlin to estimate the betas and the HRF params
      options = optimset('lsqnonlin'); % The default options for the lsqnonlin function
      betasPlusHRF = lsqnonlin(@calcModel,initialVals,lb,ub,options,roiTC,DM,TR,concatInfo);

      % given those values, calculate r2 etc and save
      glm = getGLM(betasPlusHRF,roiTC,DM,TR,concatInfo);
      saveGLMestHRF(subj,expt,model,ROI,glm) % subfunction ;
    
    else % if don't redo, then load
      glm = loadGLMestHRF(subj,expt,model,ROI);
    end % if(redo)
    
    b = glm.betas; % so can display in legend
    
    % now calculate TTA for the model
    modelTTAcorrect = getTTA(glm.model,exptTiming,1,numBins);
    modelTTAwrong = getTTA(glm.model,exptTiming,-1,numBins);
    
    % now plot
    plotNum = plotNum + 1; %goes across and then down
    figure(figNum), subplot(numRoi,numModel,plotNum);
    errorbar(-4:2:30, dataTTAcorrect(:,1),dataTTAcorrect(:,2),'r.-.');
    hold on;
    legendHandle(1) = errorbar(-4:2:30,modelTTAcorrect(:,1),modelTTAcorrect(:,2),'k.-');
    legendStr{1} = sprintf('%.2f %.2f %.2f', b(1),b(2),b(3));
    %    xlabel('time');
    xlim([-4 30]);
    %    ylabel('MRI signal');
    ylim([-0.2 0.4])
    title(sprintf('corr %s %s %0.2f',ROI,modnik,glm.r2));
    legend(legendHandle,legendStr);
    
    clear initialVals lb ub betasPlusHRF options glm modelTTAcorrect modelTTAwrong

  end
  clear roiTC dataTTAcorrect dataTTAwrong
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vals lb ub]  = setInitialVals(DM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vals(1) = 6; % rdelay;
vals(2) = 16; % udelay;

numBetas = size(DM,2);
for iB = 1:numBetas
  vals(iB+2) = 0; % initialize to 0;
end

% set lower bounds
lb(1) = 2; % rdelay;
lb(2) = 5; % udelay;
for iB = 1:numBetas
  lb(iB+2) = -2;
end

% set upper bounds
ub(1) = 12; % rdelay;
ub(2) = 26; % udelay;
for iB = 1:numBetas
  ub(iB+2) = 2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y] = calcModel(x,roiTC,DMlong,TR,concatInfo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% x are the input vals for the function
% y are the error values
% lsqnonlin will minimize the sum of squares of these errors by changing
% values in the input values, until some criteria is reached

hrfParams.rdelay = x(1);
hrfParams.udelay = x(2);
hrfParams.rdispersion = 1; % don't try to estimate;
hrfParams.udispersion = 1; % don't try to estimate;
hrfParams.tmax = 16;

hrf = spmHRF_so(TR,hrfParams);

betas(:,1) = x(3:end);

%convolve DM by HRF one run at a time and stack


% get relevant concatInfo
runTransition = concatInfo.runTransition;
hipassfilter = concatInfo.hipassfilter;

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


% calculate model given input betas and this scm from the input hrf vals
model = scm * betas;

% feed error back to lsqnonlin to minimize
y = roiTC - model; 

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function glm = getGLM(x,roiTC,DMlong,TR,concatInfo);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hrfParams.rdelay = x(1);
hrfParams.udelay = x(2);
hrfParams.rdispersion = 1; % don't try to estimate;
hrfParams.udispersion = 1; % don't try to estimate;
hrfParams.tmax = 16;

hrf = spmHRF_so(TR,hrfParams);
glm.hrf = hrf; % save the hrf

glm.betas(:,1) = x(3:end);

%convolve DM by HRF one run at a time and stack

% get relevant concatInfo
runTransition = concatInfo.runTransition;
hipassfilter = concatInfo.hipassfilter;


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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveGLMestHRF(subj, expt, model, ROI, glm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch expt
  case{'detect'}
    glmdir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/roiBetas/';
  case{'memory'}
    glmdir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/roiBetas/';
end

dataName = [glmdir 'glm_' subj '_' expt '_' ROI '_' model '_estHRF'];

save(dataName, 'glm');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function glm = loadGLMestHRF(subj, expt, model, ROI)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch expt
  case{'detect'}
    glmdir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/roiBetas/';
  case{'memory'}
    glmdir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/roiBetas/';
end

dataName = [glmdir 'glm_' subj '_' expt '_' ROI '_' model '_estHRF'];

load(dataName)
if(exist('glm'))
  glm = glm;
else
  glm = [];
end





  

