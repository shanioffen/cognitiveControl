function bootstrap = bootstrap_master(subj,expt)

% syntax: bootstrap = bootstrap_master(subj,expt)
%
% this code calls the bootstrap DM maker to get a
% set of random trial times, and then estimates the 
% betas based on those random times.
%
% it does this 10,000 per subj/expt and stores
% the results in a matrix with the following
% dimensions:
% expt x subj x ROI x beta x numBoots
% 
% since you call this once for each subj/expt,
% the output is a matrix of size
% 1 x 1 x 10 x 3 x 10,000

% 04-06-07 new ROI list

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HRFlength = 21; % better to let the model go out to end, but then for display only go to 18 seconds


% check inputs and set defaults:
if nargin < 2
    disp('Must enter a subject and experiment')
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize some values
TR = 2; % TR is 2 sec
numTpnts = 120;

redoData = 0; % whether or not to re-preprocess the data

ROIlist = {'V1V2V3_restrict','V1_restrict','V2_restrict', 'V3_restrict', 'V3AB_restrict', 'V4_restrict', 'V5_restrict', 'V7_restrict', 'LO1_restrict', 'LO2_restrict', 'VO1_restrict'}; 

betaList = {'s1', 'd', 's2'};

switch expt
    case 'memory'
	exNum = 'P5';
    case 'detect'
	exNum = 'P1';
    case 'memDetect'
	exNum = 'P3';
    case 'vertical'
	exNum = 'P4';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the output matrix
bootstrap = NaN*ones(1,1,length(betaList),length(ROIlist),1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load HRF parameters per subject

% hrf - loads bestfit. use the one for V1V2V3 combined
hrfName = ['SPM_HRF_params_' subj '_V1V2V3_restrict_' num2str(HRFlength)];
load(['/Local/Users/shani/fMRI_data/GLM_output/HRF_est/parameters/' hrfName])
fitparams = bestfit.params;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get HRF estimate using spm code; will convolve with this below

% We used the measured HRF to estimate four of the six parameters
% used by the function spm_hrf. For the other two, we'll use the 
% defaults
% The TR for the HRF expt is 1.5, but need to resample at TR of data
  
timeVec = 0:TR:30; % by 30 seconds, HRF should have returned to zero
p = [fitparams(1) fitparams(2) fitparams(3) fitparams(4) 6 0 timeVec(end)];
[hrf,p]=spm_hrf(TR,p);
hrf = hrf*(1.5/2); % need to rescale to account for the resampling because
                   % the HRF expt was at a different TR (1.5) than the data expt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate a random DM for the bootstrap
sessionDesignMatrix = mkDMbootstrap_master(subj,expt);
numRuns = size(sessionDesignMatrix,1)/numTpnts;

% convolve it with the HRF to get a model matrix:
convFiltDM = [];
for scan = 1:numRuns
  % Each scan is 120 time points, so go through the matrix
  % 120 time points at a time
  DM=sessionDesignMatrix((scan-1)*numTpnts+1:scan*numTpnts,:);

  % need to convolve each column separately
  temp = conv(DM(:,1),hrf);
  convDM(:,1) = temp(1:numTpnts); clear temp
  temp = conv(DM(:,2),hrf);
  convDM(:,2) = temp(1:numTpnts); clear temp
  temp = conv(DM(:,3),hrf);
  convDM(:,3) = temp(1:numTpnts); clear temp
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % preprocess like data:
  % Don't divide by mean, because already in % signal change
  % (since the HRF model is in % signal change)
  % filter each column (this also removes the mean so it's mean 0)
  tempDM(:,1) = FilterF([0.05 0.5],convDM(:,1));
  tempDM(:,2) = FilterF([0.05 0.5],convDM(:,2));
  tempDM(:,3) = FilterF([0.05 0.5],convDM(:,3));
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % stack to get a single DM per experiment per subject
  convFiltDM = [convFiltDM; tempDM];   
  clear tempDM     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now go through the ROIs and estimate the betas
for nROI = 1:length(ROIlist)
  ROI = ROIlist{nROI};
  filtDataName = ['filteredData_' expt '_' subj '_' ROI]; % the already-processed data
  
  %%%%%%%%%%%%%%
  % missing V5 for LM and RS, so skip those
  if and(strcmp(subj,'RS'),strcmp(ROI,'V5_restrict'))|and(strcmp(subj,'LM'),strcmp(ROI,'V5_restrict'))
    continue
  else

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Need to preprocess the data scan by scan and then stack
    % or else just load in the already processed data
    
    if(redoData)
      dataMatrix = [];
      
      % data - this loads sessionDataMatrix
      dataName = ['Raw_' expt exNum '_' subj '_' ROI '_DT_none'];
      load(['/Local/Users/shani/fMRI_data/GLM_output/' expt '/' dataName])
      
      for scan = 1:numRuns
	% Each scan is 120 time points, so go through the big scan matrices
	% 120 time points at a time
	dataM=sessionDataMatrix((scan-1)*numTpnts+1:scan*numTpnts,:);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% preprocess the data:
	% 
	% divide by mean and subtract 1:
	% This has already been done by mrLoadRet when extracting
	
	% filter to detrend
	dataM(:,1) = FilterF([0.05 0.5],dataM(:,1));
	
	% stack to get a single dataMatrix per experiment per subject
	dataMatrix = [dataMatrix; dataM];   
	clear dataM     
	
      end % for scan = 1:numRuns
      save(['/Local/Users/shani/fMRI_data/GLM_output/filteredData/' filtDataName],'dataMatrix')
    else
      load(['/Local/Users/shani/fMRI_data/GLM_output/filteredData/' filtDataName])
    end % if(redoData)
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % add a column of ones to get rid of the mean
    onesCol = ones(length(dataMatrix),1);
    convFiltDM(:,4)=onesCol;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % regress to get betas
    betas = convFiltDM\dataMatrix;
    
%    % check that the mean was near zero:
%   if betas(4)>.006
%      [num2str(betas(4)) ' mean for expt ' expt ' subj ' subj ' ROI ' ROI]
%    end
    
    % stick the betas in the output matrix and go on to the next ROI
    bootstrap(1,1,:,nROI,1) = betas(1:3); % leave out the de-mean beta
  end % end of skipping LM and RS for V5
end % end of going through ROIs - for nROI = 1:length(ROIlist)

  

    
