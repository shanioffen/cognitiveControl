function estimateAll(redo)

% the idea is to estimate the HRF for each subj/ROI using the data from both expts
% since you can't use a different HRF for different experiments:

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
  
if nargin == 0
  redo = 1; % whether to recalculate the GLMs or just load them;
end

figNum = figure; close(figNum);

roiList = {'leftDMsPCS', 'rightDMsPCS', 'leftDLsPCS', 'rightDLsPCS'};
numRoi = length(roiList);

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


for iSubj = 1:numSubj
  subj = subjList{iSubj};
  
  % load data and DM and do GLM plus estimate HRF
  if(redo)
  
    % initialize data
    for iExpt = 1:numExpt
      expt = exptList{iExpt};
      switch expt
        case{'detect'}
          datadir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/roiTseries/';
        case{'memory'}
          datadir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/roiTseries/';
      end
      
      % load model DM (for now only look at trials separated by correct/incorrect)
      DMpart{iExpt} = getDM(subj,expt,model); % subfunction;
      
      %initialize the data:
      datapart{iExpt} = NaN*ones(length(DMpart{iExpt}),numRoi); 
      % load data for all ROIs
      for iRoi = 1:numRoi
        ROI = roiList{iRoi};
        
        % load ROI time course and concatInfo
        % this loads two variables: roiTseries and concatInfo
        load([datadir subj '_' expt '_' ROI '.mat']);
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
    
    % need to stack the data for each ROI
    data = NaN*ones(length(datapart{1})+length(datapart{2}),numRoi);
    data(1:length(datapart{1}),:) = datapart{1};
    data((length(datapart{1})+1):end,:) = datapart{2};
    
    % go through ROIs again and calculate GLM and HRF
    for iRoi = 1:numRoi
      ROI = roiList{iRoi};
      roiTC = data(:,iRoi);
      
      % initialize variables
      [initialVals lb ub] = setInitialVals(DM); % subfunction
      
      % do lsqnonlin to estimate the betas and the HRF params
      options = optimset('lsqnonlin'); % The default options for the lsqnonlin function
      betasPlusHRF = lsqnonlin(@calcModel,initialVals,lb,ub,options,roiTC,DM,TR,concatInfoParts{1});
      
      % given those values, calculate r2 etc and save
      glm = getGLM(betasPlusHRF,roiTC,DM,TR,concatInfoParts);
      saveGLMestHRF(subj,model,ROI,glm) % subfunction ;
      clear glm betasPlusHRF
    end % going through ROI and doing GLM
  end % if redo

  for iRoi = 1:numRoi
    
    ROI = roiList{iRoi};
   
    % load GLM data
    glm = loadGLMestHRF(subj,model,ROI);
    r2 = glm.r2;
    hrf = glm.hrf;
    dpi = (glm.betas(8) - glm.betas(2))/(abs(glm.betas(2))+abs(glm.betas(8)));
    % that's d-mem minus d-det over the sum of absolute values
    
    for iExpt = 1:numExpt
      expt = exptList{iExpt};
      switch expt
        case{'detect'}
          datadir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/roiTseries/';
        case{'memory'}
          datadir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/roiTseries/';
      end
      
      %load data for this ROI/expt
      % this loads two variables: roiTseries and concatInfo
      load([datadir subj '_' expt '_' ROI '.mat']);
      
      % get average ROI time course
      roiTC = mean(roiTseries,2);
      
      % get Expt timing
      exptTiming = getExptTiming(subj,expt); %subfunction      

      % calculate TTA of data for plotting:
      dataTTAcorrect = getTTA(roiTC,exptTiming,1,numBins); % get TTA for correct
      dataTTAwrong = getTTA(roiTC, exptTiming,-1,numBins);
      
      % now calculate TTA for the model
      if iExpt == 1
        modelData = glm.model(1:glm.numRunsDetect*numTpnts);
        b = glm.betas(1:end/2);
      elseif iExpt == 2
        modelData = glm.model((glm.numRunsDetect*numTpnts)+1:end);
        b = glm.betas(((end/2)+1):end);
      end
      
      modelTTAcorrect = getTTA(modelData,exptTiming,1,numBins);
      modelTTAwrong = getTTA(modelData,exptTiming,-1,numBins);
      
      % get epoch-based R2
      [epochR2 resid] = calcVarAccnt(dataTTAcorrect, modelTTAcorrect);
    
      % now plot
      plotNum = iSubj + numSubj*(iExpt-1);
      figure(iRoi+figNum), subplot(numExpt+1,numSubj,plotNum);
      errorbar(-4:2:30, dataTTAcorrect(:,1),dataTTAcorrect(:,2),'r.-.');
      hold on;
      legendHandle(1) = errorbar(-4:2:30,modelTTAcorrect(:,1),modelTTAcorrect(:,2),'k.-');
      legendStr{1} = sprintf('%.2f %.2f %.2f', b(1),b(2),b(3));
      %    xlabel('time');
      xlim([-4 30]);
      %    ylabel('MRI signal');
      ylim([-0.2 0.4])
      title(sprintf('%s dpi = %0.2f r2 = %.2f',expt,dpi,epochR2));
      legend(legendHandle,legendStr);
      
      clear roiTseries roiTC exptTiming dataTTAcorrect dataTTAwrong modelData b ;
      clear modelTTAcorrect modelTTAwrong legendHandle legendStr;
      clear epochR2 resid;
    end % going through expt
    % plot the HRF being used
    plotNum = iSubj + 2*numSubj;
    figure(iRoi + figNum), subplot(numExpt+1,numSubj,plotNum);
    lHand = plot(hrf);
    lStr = ROI;
    title(sprintf('HRF est for %s, r2 = %.2f',subj,r2));
    legend(lHand,lStr);
    clear glm hrf r2
  end
  clear roiTC dataTTAcorrect dataTTAwrong dpi ROI
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
hrfParams.tmax = 16;

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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveGLMestHRF(subj, model, ROI, glm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
glmdir = '/Users/shani/NYU/fMRI_data/GLM_output/RoB_ROIanalysis/';
dataName = [glmdir 'glm_' subj '_' ROI '_' model '_estHRF_bothExpt'];

save(dataName, 'glm');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function glm = loadGLMestHRF(subj, model, ROI)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
glmdir = '/Users/shani/NYU/fMRI_data/GLM_output/RoB_ROIanalysis/';


dataName = [glmdir 'glm_' subj '_' ROI '_' model '_estHRF_bothExpt'];

load(dataName)
if(exist('glm'))
  glm = glm;
else
  glm = [];
end





  

