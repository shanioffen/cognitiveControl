function allparams =  estimateHRFs(restrict)
% function allparams =  estimateHRFs(restrict)
% the idea is to estimate the HRF for each subj/ROI using the data from both expts
% since you can't use a different HRF for different experiments:
% if enter 1 for restrict (default), use data from restricted ROIs
% where restriction is based on r2 model fit using canonical hrf
% the HRF parameters are saved out in a variable glm.hrfParams
% in the folder glmdir = '/Users/shani/NYU/fMRI_data/GLM_output/RoB_ROIanalysis/';
% and the dataName = [glmdir 'glm_' subj '_' ROI '_' model '_estHRF_bothExpt'];


% The main goal of this code is to calcualate estimated HRF parameters
% that can then be used just like the other HRF parameters to estimate models,
% but I also save out the rest of the glm

% This code concatenates the data for both experiments, to estimate
% a single HRF for both, but separate betas for each.
% stick with one model = S2plus (models S2 with 2 TRs)
% if just want to plot without calculating, give argument of 0
% (eg estimateAll(0))
  
if nargin ==0
  restrict = 1;
end

figNum = figure; close(figNum);

ROIlist = {'leftAntCS','leftDLsPCS','leftDMsPCS','leftIPCS','leftPostIPS',...
           'rightAntCS','rightDLsPCS','rightDMsPCS','rightIPCS','rightPostIPS'};
numRoi = length(ROIlist);

exptList = {'detect','memory'};
numExpt = length(exptList);

subjList = {'JG','DS','LM','RS','SO'};
numSubj = length(subjList);

% fixing one model but leave all this here just in case...
modelList = {'correctTrial','correctTrialS2plus','s2plusLong','correctTrialS2trip','S2tripLong'};  
modelNickname = {'S2single','S2plus','S2plusLong','S2trip','S2tripLong'};
numModel = length(modelList);
model = modelList{2}; modelnik = modelNickname{2};

TR = 2;
numBins = 1; % not plotting all durations;
numTpnts = 120;

glmDir = '/Users/shani/NYU/fMRI_data/GLM_output/RoB_ROIanalysis/';

for iSubj = 1:numSubj
  subj = subjList{iSubj};
  
  % initialize data
  for iExpt = 1:numExpt
    expt = exptList{iExpt};
    switch expt
      case{'detect'}
        datadir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/roiTseries/';
      case{'memory'}
        datadir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/roiTseries/';
    end
    
    % load DM (for now only look at trials separated by correct/incorrect)
    DMpart{iExpt} = getDM(subj,expt,model); % subfunction;
    
    %initialize the data:
    datapart{iExpt} = NaN*ones(length(DMpart{iExpt}),numRoi); 
    % load data for all ROIs
    for iRoi = 1:numRoi
      ROI = ROIlist{iRoi};
      
      % load ROI time course and concatInfo
      % this loads two variables: roiTseries and concatInfo
      dataFileName = [datadir subj '_' expt '_' ROI];
      if restrict
        dataFileName = [dataFileName '_restCan25.mat'];
      else
        dataFileName = [dataFileName '.mat'];
      end
      
      load(dataFileName)
      concatInfoParts{iExpt} = concatInfo;
      
      % get average ROI time course and store by ROI
      roiTC = mean(roiTseries,2);
      datapart{iExpt}(:,iRoi) = roiTC;
      clear roiTC roiTseries ROI;
    end %going through ROI
  end % going through experiments to get DM and data
    
  % need to create a DM that has upper left hand corner the DM for detect,
  % and lower right hand corner the DM for memory
  DM = zeros(size(DMpart{1},1)+size(DMpart{2},1),2*size(DMpart{1},2));
  DM(1:size(DMpart{1},1),1:size(DMpart{1},2)) = DMpart{1};
  DM(size(DMpart{1},1)+1:end,size(DMpart{1},2)+1:end) = DMpart{2};
  clear DMpart;
  
  % average the data across ROIs and concatenate
  avgeData{1} = mean(datapart{1},2);
  avgeData{2} = mean(datapart{2},2);
  clear datapart;
  
  data = [avgeData{1}; avgeData{2}];
  %check I did this right
  if length(data) ~= length(avgeData{1}) + length(avgeData{2})
    disp('didnt concatenate data correctly');
    keyboard;
  end
  
  % use lsqnonlin to estimate HRF
  % initialize variables
  [initialVals lb ub] = setInitialVals(DM); % subfunction
  
  % do lsqnonlin to estimate the betas and the HRF params
  options = optimset('lsqnonlin'); % The default options for the lsqnonlin function;
  disp(sprintf('starting lsqnonlin, subj %s',subj))
  betasPlusHRF = lsqnonlin(@calcModel,initialVals,lb,ub,options,data,DM,TR,concatInfoParts{1});
  disp('done')
  % save the HRF parameters
  estParams.rdelay = betasPlusHRF(1);
  estParams.udelay = betasPlusHRF(2);
  
  % also return them in a matrix
  allparams(1:2,iSubj) = betasPlusHRF(1:2);
  
  dataName = [glmDir 'HRFest_' subj];
  if restrict
    dataName = [dataName '_restCan25'];
  end
  
  save(dataName,'estParams');
  clear hrfParams rdelay udelay
  
  % given those values, calculate r2 etc to plot
  glm = getGLM(betasPlusHRF,data,DM,TR,concatInfoParts);
  r2 = glm.r2;
  hrf = glm.hrf;

  % now plot to check it out  

  detectTiming = getExptTiming(subj,'detect'); 
  memoryTiming = getExptTiming(subj,'memory');
  
  % calculate TTA of data for plotting:
  dataTTAdetect = getTTA(avgeData{1}, detectTiming,1,1); 
  dataTTAmemory = getTTA(avgeData{2}, memoryTiming,1,1);
  
  modelDetect = glm.model(1:glm.numRunsDetect*numTpnts);
  modelMemory = glm.model((glm.numRunsDetect*numTpnts)+1:end);
  
  modelTTAdetect = getTTA(modelDetect,detectTiming,1,1);
  modelTTAmemory = getTTA(modelMemory,memoryTiming,1,1);
  
  % get epoch-based R2
  [epochR2d resid] = calcVarAccnt(dataTTAdetect, modelTTAdetect);
  [epochR2m resid] = calcVarAccnt(dataTTAmemory, modelTTAmemory);    
  
  % now plot
  figure(figNum), subplot(3,5,iSubj)
  errorbar(-4:2:30, dataTTAdetect(:,1),dataTTAdetect(:,2),'r.-');
  hold on;
  legendHandle(1) = errorbar(-4:2:30,modelTTAdetect(:,1),modelTTAdetect(:,2),'k.-');
  legendStr{1} = ('DETECT');
  xlim([-4 30]);
  ylim([-0.2 0.4])
  title(sprintf('%s r2 = %.2f',subj,epochR2d));
  legend(legendHandle,legendStr);
  clear legendHandle legendStr
  
  figure(figNum), subplot(3,5,iSubj+numSubj)
  errorbar(-4:2:30, dataTTAmemory(:,1),dataTTAmemory(:,2),'r.-');
  hold on;
  legendHandle(1) = errorbar(-4:2:30,modelTTAmemory(:,1),modelTTAmemory(:,2),'k.-');
  legendStr{1} = ('MEMORY');
  xlim([-4 30]);
  ylim([-0.2 0.4])
  title(sprintf('%s r2 = %.2f',subj,epochR2m));
  legend(legendHandle,legendStr);
  clear legendHangle legendStr
  
  % plot the HRF being used
  figure(figNum), subplot(3,5,iSubj+2*numSubj)
  lh(1) = plot(hrf,'k');
  ls{1} = 'est';
  hold on;
  lh(2) = plot(spm_hrf(2),'b');
  ls{2} = 'can';
  title(sprintf('r2 = %.2f',r2));
  legend(lh,ls);
  clear lh ls;
  
  clear dataTTAdetect dataTTAmemory modelTTAdetect modelTTAmemory;
  clear epochR2d epochR2m detectTiming memoryTiming;
  clear avgeData modelDetect modelMemory glm r2 hrf


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
hrfParams.tmax = 32;

hrf = spmHRF_so(TR,hrfParams);

betas(:,1) = x(3:end);

%convolve DM by HRF one run at a time and stack


% get relevant concatInfo
hipassfilter = concatInfo.hipassfilter{1}; % same for all runs and all expts so just take first;

numTpnts = 120;
numRuns = length(DMlong)/numTpnts;

scm = [];
for runnum = 1:numRuns
  DM = [];
  % convolve DM with HRF
  DM = DMlong(((numTpnts*(runnum-1))+1):(numTpnts*runnum),:);
  m = convn(DM, hrf(:,1));
  m = m(1:length(DM),:);
  % remove mean 
  m = m-repmat(mean(m), size(m,1), 1);
  % apply the same filter as original data
  m = real(ifft(fft(m) .* repmat(hipassfilter', 1, size(m,2)) ));
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
hrfParams.tmax = 32;

hrf = spmHRF_so(TR,hrfParams);
glm.hrf = hrf; % save the hrf;
glm.hrfParams = hrfParams; % save the params so can be used again

glm.betas(:,1) = x(3:end);

%convolve DM by HRF one run at a time and stack

% get relevant concatInfo
hipassfilter = concatInfo{1}.hipassfilter{1}; % same for all runs and all expts so just take first;

% keep track of how many runs for detection vs memory
glm.numRunsDetect = concatInfo{1}.n;
glm.numRunsMemory = concatInfo{2}.n;
% make sure they add up:
numTpnts = 120;
numRuns = length(DMlong)/numTpnts;

if numRuns ~= glm.numRunsDetect + glm.numRunsMemory
  disp('estimateAll.m:getGLM: the number of runs dont add up')
  return
end

%convolve DM by HRF one run at a time and stack
scm = [];
% test how long this takes
for runnum = 1:numRuns
  DM = [];
  % convolve DM with HRF
  DM = DMlong(((numTpnts*(runnum-1))+1):(numTpnts*runnum),:);
  m = convn(DM, hrf(:,1));
  m = m(1:length(DM),:);
  % remove mean 
  m = m-repmat(mean(m), size(m,1), 1);
  % apply the same filter as original data
  m = real(ifft(fft(m) .* repmat(hipassfilter', 1, size(m,2)) ));
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








  

