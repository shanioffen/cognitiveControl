function [betas] = estimateBetas(subj,expt,HRFnum)
% function estimateBetas(subj,expt,[HRFnum])
% calculate betas for given Subj,Expt for all ROIs
% plots time courses
% Model chosen is S2plus (modeling the end of hte trial with 2TRs)
%
% can use canonical HRF or estimated HRF.
% for canonical, enter HRFnum==1 or leave blank;
% for estimated, enter HRFnum==2.
% to use estimated HRF, need to first estimate it using code 'estimateHRFs.m'
  
redo = 1; % whether to recalculate the GLMs or just load them and show graphs;
doplot = 1; % whether to graph the timecourses or just calculate the betas

if nargin==2
  HRFnum = 1;
end

hrfList = {'Can','Est','Vis'};
whichHRF = hrfList{HRFnum};

ROIlist = {'leftAntCS','leftDLsPCS','leftDMsPCS','leftIPCS','leftPostIPS',...
           'rightAntCS','rightDLsPCS','rightDMsPCS','rightIPCS','rightPostIPS'};
numRoi = length(ROItags);

model = {'correctTrialS2plus'};
modelnik = {'S2plus'};

switch expt
  case{'detect'}
    datadir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/roiTseries/';
    % if doing fixed effects, need to keep track of how many runs each sub did, for HRF 
    % *** order of scans is JG, DS, LM, RS, SO ********    
    subjRunList = [10 17 10 10 20];
  case{'memory'}
    datadir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/roiTseries/';
    % if doing fixed effects, need to keep track of how many runs each sub did, for HRF 
    % *** order of scans is JG, DS, LM, RS, SO ********
    subjRunList = [11 14 12 11 10];
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
  [runTransition hipassfilter] = getConcatInfo(subj,expt);
  
  % calculate TTA for plotting:
  dataTTAcorrect = getTTA(roiTC,exptTiming,1); % get TTA for correct
  dataTTAwrong = getTTA(roiTC, exptTiming,-1);
  
  % do GLM on that time course if haven't already
  if(redo)
    % load model DM (for now only look at trials separated by correct/incorrect)
    DMlong = getDM(subj,expt,model); % subfunction
    
    % load the HRF that will be used to convolve the DM
    hrfParams = loadHRFparams(subj,whichHRF);
    hrfParams.tmax = 30;
    hrf = spmHRF_so(2,hrfParams); % separate fxn;
    glm.hrf = hrf; % save it to structure
    glm.hrfParams = hrfParams; %keep the paramters used to make the HRF
    
    %convolve DM by HRF one run at a time and stack
    allscm = [];
    
    for runnum = 1:size(runTransition,1)
      % determine which subjects runs these are
      if length(subjRunList) % if there are more than one subject
        for icount = length(subjRunList):-1:1
          if runnum <= sum(subjRunList(1:icount))
            whichSub = icount;
          end
        end
      else
        whichSub = 1;
      end
      
      % default values
      scm = [];
      DM = [];
      % convolve DM with HRF
      DM = d.DM(runTransition(runnum,1):runTransition(runnum,2),:);
      m = convn(DM, hrf(:,whichSub));
      m = m(1:length(DM),:);
      % remove mean 
      m = m-repmat(mean(m), size(m,1), 1);
      % apply the same filter as original data
      if ~isempty(hipassfilter)
        m = real(ifft(fft(m) .* repmat(hipassfilter{runnum}', 1, size(m,2)) ));
      end
      scm = m; % no longer doing each column separately
               % stack this run's stimcmatrix on to the last one
      allscm = [allscm;scm];
    end % going through runs
    
    glm.scm = allscm; % save the convolved DM;
    clear scm %make sure you don't call it by mistake instead of allscm
    
    % calculate the GLM and get betas
    % precalculate the normal equation (this dramatically speeds up things)
    precalcmatrix = ((allscm'*allscm)^-1)*allscm';
    % if this don't work then do pinv
    if sum(isnan(precalcmatrix(:))) == length(precalcmatrix(:))
      disp(sprintf('(estimateBetas) Using pseudo inverse to invert convolution matrix'));
      precalcmatrix = pinv(allscm);
    end
    
    glm.betas = precalcmatrix*roiTC;
    % calculate error bars, first get sum-of-squares of residual
    % (in percent signal change)
    glm.model = allscm * glm.betas;
    sumOfSquaresResidual = sum((roiTC-glm.model).^2);
    % now calculate the sum-of-squares of that error
    % and divide by the degrees of freedom (n-k where n
    % is the number of timepoints in the scan and k is 
    % the number of timepoints in all the estimated hdr)
    S2 = sumOfSquaresResidual/(length(roiTC)-size(allscm,2));
    % now distribute that error to each one of the points
    % in the hemodynamic response according to the inverse
    % of the covariance of the stimulus convolution matrix.
    glm.betaError = sqrt(diag(pinv(allscm'*allscm))*S2);
    % calculate variance accounted for by the estimated hdr
    glm.r2 = (1-sumOfSquaresResidual./sum(roiTC.^2));
    
    % save all that info
    saveGLM(subj,expt,model,ROI,glm,whichHRF) % subfunction ;
  else % if don't redo, then load
    glm = loadGLM(subj,expt,model,ROI,whichHRF);
  end % if(redo)
  
  betas = glm.betas;
  
  % now calculate TTA for the model
  modelTTAcorrect = getTTA(glm.model,exptTiming,1);
  modelTTAwrong = getTTA(glm.model,exptTiming,-1);
  
  % get epoch-based R2
  [epochR2 resid] = calcVarAccnt(dataTTAcorrect, modelTTAcorrect);
  
  % now plot
  if doplot
    plotNum = plotNum + 1; % goes across, then down
    figure(figNum), subplot(2,5,plotNum);
    legendHandle(1) = errorbar(-4:2:30, dataTTAcorrect(:,1),dataTTAcorrect(:,2),'r.-.');
    legendStr{1} = sprintf('%s %s',subj,expt);
    hold on;
    legendHandle(2) = errorbar(-4:2:30,dataTTAwrong(:,1),dataTTAwrong(:,2),'b-');
    legendStr{2} = 'Incorrect';
    
    legendHandle(2) = errorbar(-4:2:30,modelTTAcorrect(:,1),modelTTAcorrect(:,2),'k.-');    
    legendStr{3} = sprintf('r2=%0.2f, d=%0.2',glm.r2,betas(2));
    xlim([-4 30]);
    ylim([-0.2 0.4])
    title(sprintf('%s %0.2f',ROI,epochR2));
    legend(legendHandle,legendStr);
    
  end % if doplot
end %going through rois



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [runTransition hipassfilter] = getConcatInfo(subj,expt,concatInfo);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch subj
  case 'allSubj'
    runTransition = [];
    hipassfilter = concatInfo{1}.hipassfilter; % same for all
    for iSubj = 1:length(concatInfo)
      runTemp = concatInfo{iSubj}.runTransition;
      runTransition = [runTransition; runTemp];
    end
  otherwise
    runTransition = concatInfo.runTransition;
    hipassfilter = concatInfo.hipassfilter;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hrfParams = loadHRFparams(subj,whichHRF)
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
        dataName = [glmHRFdir 'glm_' subj '_estHRF_RoB'];      
        load(dataName); % this loads variable glm with field glm.hrfParams
        hrfParams.rdelay(iSub,1) = glm.hrfParams.rdelay;
        hrfParams.udelay(iSub,1) = glm.hrfParams.udelay;
        hrfParams.udispersion(iSub,1) = 1;
        hrfParams.rdispersion(iSub,1) = 1;
        clear glm;
      end
    else % just load for individual subject if not doing fixed effects
      dataName = [glmHRFdir 'glm_' subj '_estHRF_RoB'];      
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

