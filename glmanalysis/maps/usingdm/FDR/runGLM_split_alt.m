% 2007August19
% changed to have new HRF
% 2008March10 - expanded to use different models (Added the 'which')
% 2008March11 - having a problem running in loops so running by subj/expt instead, see if that works
% 2008March17 - added ability to blur when concatenating
% 2008March19 - added ability to run all subjects for fixed effects (needed to fix how HRFs are loaded)
% 2008 May 5 - added ability to use HRF as estimated in frontal cortex using estimateHRFall.m in ROIanalysis
% 2008 May 17 - added FDR overlay


function runGLM_split_alt(subj,expt)
% whichModel = {'default','corrDelay','corrTrial','split','corrDelay_S2plus','corrTrial_S2plus'};
% (the difference between models:) 
% default uses 3 predictors: s1, delay, s2 (treats all trials the same)
% correctDelay uses 4 predictors: s1, delay_correct, delay_wrong, s2
% correctTrial uses 6 predictors: s1_correct, delay_correct, s2_correct, and s1,d,s2 _wrong
% s2plus means model S2 with 2 TRs instead of 1 TR
% split is like correct Trial but analyzes half the data at a time to test for robustness

hrfList = {'Vis','Can','Est'};
numHrf = length(hrfList);
whichHRF = hrfList{3};

clear view
clear global MLR
clear global mrDEFAULTS

baseDir = setBaseDir(expt,subj);
dataDir = [baseDir subj '_' expt];
cd(dataDir)
mrGlobals;
mrSetPref('overwritePolicy','Merge');
mrSetPref('verbose','No');
view = newView('Volume');

modelList = {'default','correctDelay','correctTrial','correctDelayS2plus','correctTrialS2plus','correctTrialS2minus'};
modelNickname = {'all','byDelay','byTrial','byDelayS2plus','byTrialS2plus','byTrialS2minus'};
modelNum = 5;

for iModel = modelNum
  model = modelList{iModel};
  modelnik = modelNickname{iModel};
  
  % set parameters to run GLM
  groupNames = viewGet(view,'groupNames');
  params.saveName = ['FDRglm_' modelnik 'hrf_' whichHRF];
  params.hrfModel = 'spmHRF_so';
  params.trSupersampling = 1;
  params.groupName = 'ConcatBlur'; 
  params.whichModel = model;
  
  % now explicitly setting the DM
  DM = getDM(subj,expt,model);
  params.DM = DM;
  
  if isempty(params.DM{1}),mrWarnDlg('(runGLM_split) No DM found');keyboard,end
  
  params.subjRunList = [];
  
  paramsInfo = {...
      {'groupName',groupNames,'Name of group from which to do eventRelated analysis'},...
      {'saveName',params.saveName,'File name to try to save as'},...
      {'hrfModel',params.hrfModel,'Name of the function that defines the hrf used in glm'},...
      {'trSupersampling', params.trSupersampling, 'minmax=[1 100]', 'TR supersampling factor (1=no supersampling) reulting design matrix will be downsampled afterwards'}};
  
  % load the hrf parameters for the subject (subfunction)
  params.hrfParams = loadHRFparams(subj,whichHRF);
  params.hrfParams.description = params.hrfModel;
  params.hrfParams.incDeriv = 0;
  params.hrfParams.paramInfo  = feval(params.hrfModel, 'params');
  
  view = viewSet(view,'groupName',params.groupName);
  groupNum = viewGet(view,'currentGroup');
  
  params.scanNum = inputScanNum;
  useDefault = 1;
  params.scanParams = getEventRelatedParams(view,params,useDefault);
  for i = 1:length(params.scanNum)
    params.scanParams{params.scanNum(i)}.scanNum = params.scanNum(i);
    params.scanParams{params.scanNum(i)}.description = ['GLM analysis of Concatenation: ' num2str(i)];
    params.scanParams{params.scanNum(i)}.hdrlen = 30; % how far out to go, in seconds
  end

  disp('(runGLM_split) calling mrLR_FDRoverlay for all runs')
  [view d] = mrLR_FDRoverlay(view,params);

end % iModel


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function baseDir = setBaseDir(expt,subj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set directories
if strcmp(subj,'allSubj')
  switch expt
    case 'memory'
      baseDir = '/Volumes/extraDrive/NYU/fMRI_data/Memory/eventRelated/';
    case 'detect'
      baseDir = '/Volumes/extraDrive/NYU/fMRI_data/Detection/detectionEvent/';
    case 'memDetect'
      baseDir = '/Volumes/extraDrive/NYU/fMRI_data/MemDetect/';
    case 'vertical'
      baseDir = '/Volumes/extraDrive/NYU/fMRI_data/Vertical/';
  end % switch expt
else
  switch expt
    case 'memory'
      baseDir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/';
    case 'detect'
      baseDir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/';
    case 'memDetect'
      baseDir = '/Volumes/extraDrive/NYU/fMRI_data/MemDetect/';
    case 'vertical'
      baseDir = '/Volumes/extraDrive/NYU/fMRI_data/Vertical/';
  end % switch expt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = stim2mrLR(v,model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

groupNum = viewGet(v,'groupNum','MotionComp');
nScans = viewGet(v,'nScans',groupNum);
for iScan=2:nScans-1 % first and last scans are localizers
  stimFilename = ['stimFiles/stimFile_' model '_' num2str(iScan-1) '.mat'];
  v = viewSet(v,'stimFilename',stimFilename,iScan,groupNum);
end

saveSession;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function view = junkAndCat(view,blur, split)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% change junk frames and concatenate

global MLR
for groupNum = [viewGet(view,'groupNum','Raw') viewGet(view,'groupNum','MotionComp')]; % changing for raw and motion corrected
  for iScan = 1:viewGet(view,'nScans',groupNum)
    MLR.groups(groupNum).scanParams(iScan).junkFrames = 5; %0; % 5; 
    MLR.groups(groupNum).scanParams(iScan).nFrames =  120; %125; % 120; 
    MLR.groups(groupNum).scanParams(iScan).totalFrames = 125;
  end
  saveSession;
end

% get params
% [view, params] = concatTSeries(view,[],'justGetParams=1','defaultParams=1');

% set params
params.groupName = 'MotionComp';
if blur, params.newGroupName = 'ConcatBlur'; else,  params.newGroupName = 'Concatenation'; end
params.description = 'Concatenation of [x...x] with hipass filtering'; 
params.filterType = 1; 
params.filterCutoff = 0.0100;
params.percentSignal = 1;
params.warp = 1; params.warpBaseScan = 2; 
params.warpInterpMethod = 'linear';
params.blur = blur;

if split == 0 % concatenate all the data
  params.scanList = 2:(viewGet(view,'nScans',viewGet(view,'groupNum',params.groupName))-1);
  view = concatNblur(view,params);
else % split the data in half, every other run
  for iSplit = 1:2
    params.scanList = 1+iSplit:2:(viewGet(view,'nScans',viewGet(view,'groupNum',params.groupName))-1);
    view = concatNblur(view,params);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hrfParams = loadHRFparams(subj,whichHRF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hrfDir = '/Users/shani/NYU/fMRI_data/GLM_output/HRF_est/parameters/';

% estimated HRF using detect and memory both, since need to use
% same HRF for both experiments.
glmDir = '/Users/shani/NYU/fMRI_data/GLM_output/RoB_ROIanalysis/';

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
        dataName = [glmDir 'HRFest_' sub];
        load(dataName); % this loads variable hrfParams;
        hrfParams.rdelay(iSub,1) = estParams.rdelay;
        hrfParams.udelay(iSub,1) = estParams.udelay;
        hrfParams.udispersion(iSub,1) = 1;
        hrfParams.rdispersion(iSub,1) = 1;
        clear glm;
      end
    else % just load for individual subject if not doing fixed effects
      dataName = [glmDir 'HRFest_' subj];
      load(dataName); % this loads variable hrfParams
        hrfParams.rdelay = estParams.rdelay;
        hrfParams.udelay = estParams.udelay;
        hrfParams.udispersion = 1;
        hrfParams.rdispersion = 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DM  = getDM(subj, expt, model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dmDir = '/Users/shani/NYU/fMRI_data/GLM_output/DM/'; %where to save the design matrices
dataName1 = [dmDir 'DM_SPLITcorrectDelayS2plus_' expt '_' subj]; % this is really corretTrial, jsut named wrong by typo in DM making code!
dataName2 = [dmDir 'DM_correctTrialS2plus_' expt '_' subj];
load(dataName1);
load(dataName2);
DMall = DM; clear DM;

if(exist('DModd'))
  DM{1} = DMall;
  DM{2} = DModd;
  DM{3} = DMeven;
else
  DM{1} = [];
  DM{2} = [];
  DM{3} = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function scanParams = getEventRelatedParams(view,params,useDefault)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TAKEN FROM eventRelatedGlmGUI.m

% make the output as long as the number of scans
scanParams = cell(1,viewGet(view,'nScans',viewGet(view,'groupNum',params.groupName)));

% check for stimfile, and if it is mgl/type then ask the
% user which variable they want to do the anlysis on
for scanNum = 1:length(params.scanNum)
  % get scan and default description
  scanInfo = sprintf('%i: %s',params.scanNum(scanNum),viewGet(view,'description',params.scanNum(scanNum)));
  description = sprintf('Event related analysis of %s: %i',params.groupName,params.scanNum(scanNum));
  % standard parameters to set
  taskVarParams = {...
      {'scan',scanInfo,'type=statictext','Description of scan to set parameters for (not editable)'},...
      {'description',description,'Event related analysis of [x...x]','Description of the analysis'}...
      {'hdrlen',25,'Length of response in seconds to calculate'}...
      {'preprocess','','String of extra commands for preprocessing. Normally you will not need to set anything here, but this allows you to do corrections to the stimvols that are calculated so that you can modify the analysis. (see wiki for details)'}...
		  };

    
  % give the option to use the same variable for all
  if (scanNum == 1) && (length(params.scanNum)>1)
    taskVarParams{end+1} = {'sameForAll',1,'type=checkbox','Use the same variable name for all analyses'};
  end
  %%%%%%%%%%%%%%%%%%%%%%%
  % now we have all the dialog information, ask the user to set parameters
  if useDefault
    scanParams{params.scanNum(scanNum)} = mrParamsDefault(taskVarParams);
  else
    scanParams{params.scanNum(scanNum)} = mrParamsDialog(taskVarParams);
  end
  % user hit cancel
  if isempty(scanParams{params.scanNum(scanNum)})
    scanParams = [];
    return
  end
  %%%%%%%%%%%%%%%%%%%%%%%
    
  % check if the varname is a cell array, then convert to a cell array
  % instead of a string this is so that the user can specify a variable
  % name like {{'varname'}}
  if (isfield(scanParams{params.scanNum(scanNum)},'varname') &&...
      isstr(scanParams{params.scanNum(scanNum)}.varname) && ...
      (length(scanParams{params.scanNum(scanNum)}.varname) > 1) && ...
      (scanParams{params.scanNum(scanNum)}.varname(1) == '{'))
    scanParams{params.scanNum(scanNum)}.varname = eval(scanParams{params.scanNum(scanNum)}.varname);
  end

  % if sameForAll is set, copy all parameters into all scans and break out of loop
  if isfield(scanParams{params.scanNum(scanNum)},'sameForAll') && ...
	scanParams{params.scanNum(scanNum)}.sameForAll
    for i = 2:length(params.scanNum)
      % set the other scans params to the same as this one
      scanParams{params.scanNum(i)} = scanParams{params.scanNum(1)};
      % change the description field appropriately for this scan num
      description = scanParams{params.scanNum(1)}.description;
      groupNameLoc = strfind(description,params.groupName);
      if ~isempty(groupNameLoc)
	description = sprintf('%s%s: %i',description(1:groupNameLoc(1)),params.groupName,params.scanNum(i));
      end
      scanParams{params.scanNum(i)}.description = description;
    end
    break
  end
  taskVarParams = {};
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function hrf = spmHRF_so(TR,params) %

% **** now this is done in a separate function, not a subfunction, so doesn't get confusing since call it from many different functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




