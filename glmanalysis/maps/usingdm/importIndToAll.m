% 2008March31 - adapted runGLM to do the all subj 4 at a time leaving out 1

function importIndToAll(expt)
% function runGLM_allSubj_by4s(expt,inputModelNum)
% whichModel = {'default','corrDelay','corrTrial','split','by4s'};
% inputModelNum = 4 for splitting, = 5 for leaving one out at a time

% whichModel = {'default','corrDelay','corrTrial','split','by4s'};
% (the difference between models:) 
% default uses 3 predictors: s1, delay, s2 (treats all trials the same)
% correctDelay uses 4 predictors: s1, delay_correct, delay_wrong, s2
% correctTrial uses 6 predictors: s1_correct, delay_correct, s2_correct, and s1,d,s2 _wrong
% split is like correct Trial but analyzes half the data at a time to test for robustness

subj = 'allSubj';
whichModel = {'correctTrial','correctTrialS2plus','correctTrialS2minus'};
blur = 1; % set to 1 in order to spatially blur the data, to 0 to leave unblurred
split = 0; % do that in individual folders

% set the directory and start mrLR
clear view
clear global MLR
clear global mrDEFAULTS
baseDir = setBaseDir(expt,subj);
dataDir = [baseDir subj '_' expt];
cd(dataDir)
mrGlobals;
mrSetPref('overwritePolicy','merge');
mrSetPref('verbose','No');

view = newView('Volume');
view = stim2mrLR(view, 'corrTrial');

% if doing fixed effects, need to keep track of how many runs each sub did, 
% so you know which HRF to apply when running the GLM
% *** order of subjects is JG, DS, LM, RS, SO ********
subjList = {'JG','DS','LM','RS','SO'};
switch expt
  case 'detect'
    subjRunList = [10 17 10 10 20];
  case 'memory'
    subjRunList = [11 14 12 11 10];
end % switch expt
  

for iSub = 1:5 % go through and get each; 
  individual = subjList{iSub};

  if 0
    view = junkAndCat(view, blur, split, iSub, subjRunList);
  end

  for iModel = 3
    model = whichModel{iModel};
    % set parameters to run GLM
    groupNames = viewGet(view,'groupNames');
    params.saveName = ['glmbyDM_' model];
    params.hrfModel = 'spmHRF_so';
    params.trSupersampling = 1;
    params.groupName = 'Individuals';
    params.whichModel = model;
    params.subjRunList = [];
    % now explicitly setting the DM
    params.DM = getDM(individual,expt,model);
    if isempty(params.DM),mrWarnDlg('(runGLMbyDM) No DM found');keyboard,end
    
    paramsInfo = {...
        {'groupName',groupNames,'Name of group from which to do eventRelated analysis'},...
        {'saveName',params.saveName,'File name to try to save as'},...
        {'hrfModel',params.hrfModel,'Name of the function that defines the hrf used in glm'},...
        {'trSupersampling', params.trSupersampling, 'minmax=[1 100]', 'TR supersampling factor (1=no supersampling) reulting design matrix will be downsampled afterwards'}};
    
    % load the hrf parameters for the subject (subfunction)
    params.hrfParams = loadHRFparams(individual);
    params.hrfParams.description = params.hrfModel;
    params.hrfParams.incDeriv = 0;
    params.hrfParams.paramInfo  = feval(params.hrfModel, 'params');
    
    view = viewSet(view,'groupName',params.groupName);
    groupNum = viewGet(view,'currentGroup');
    params.scanNum = iSub;
    
    useDefault = 1;
    params.scanParams = getEventRelatedParams(view,params,useDefault);
    
    for i = 1:length(params.scanNum)
      params.scanParams{params.scanNum(i)}.scanNum = params.scanNum(i);
      params.scanParams{params.scanNum(i)}.description = ['GLM analysis of ConcatBlur: ' num2str(i)];
      params.scanParams{params.scanNum(i)}.hdrlen = 30; % how far out to go, in seconds
    end
    
    % run the GLM
    disp('(runGLMbyDM) calling mrLR_GLM_so')
    [view d] = mrLR_GLM_so(view,params);
  
  end % going through the models
  
end % going through the subjects
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function baseDir = setBaseDir(expt,subj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set directories
if strcmp(subj,'allSubj')
  switch expt
    case 'memory'
      baseDir = '/Volumes/LaCie/NYU/fMRI_data/Memory/eventRelated/';
    case 'detect'
      baseDir = '/Volumes/LaCie/NYU/fMRI_data/Detection/detectionEvent/';
    case 'memDetect'
      baseDir = '/Volumes/LaCie/NYU/fMRI_data/MemDetect/';
    case 'vertical'
      baseDir = '/Volumes/LaCie/NYU/fMRI_data/Vertical/';
  end % switch expt
else
  switch expt
    case 'memory'
      baseDir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/';
    case 'detect'
      baseDir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/';
    case 'memDetect'
      baseDir = '/Volumes/LaCie/NYU/fMRI_data/MemDetect/';
    case 'vertical'
      baseDir = '/Volumes/LaCie/NYU/fMRI_data/Vertical/';
  end % switch expt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = stim2mrLR(v, model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

groupNum = viewGet(v,'groupNum','MotionComp');
nScans = viewGet(v,'nScans',groupNum);
for iScan=2:nScans-1 % first and last scans are localizers
  stimFilename = ['stimFiles/stimFile_' model '_' num2str(iScan-1) '.mat'];
  v = viewSet(v,'stimFilename',stimFilename,iScan,groupNum);
end

saveSession;
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function view = junkAndCat(view,blur, split, iSub, subjRunList)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% change junk frames and concatenate, leaving out one subj at a time

global MLR
% set params
params.groupName = 'MotionComp';
params.newGroupName = 'Individuals'; 
params.description = 'Concatenation of [x...x] with hipass filtering'; 
params.filterType = 1; 
params.filterCutoff = 0.0100;
params.percentSignal = 1;
params.warp = 1; 
params.warpInterpMethod = 'linear';
params.blur = blur;

if split == 0 % concatenate one subj at a time
  if iSub>1
    params.scanList = sum(subjRunList(1:iSub-1))+2:sum(subjRunList(1:iSub))+1;
  else
    params.scanList = 2:subjRunList(1)+1;
  end
  params.warpBaseScan = params.scanList(1); % choose first scan in list as base
  view = concatNblur(view,params);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hrfParams = loadHRFparams(subj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hrfDir = '/Users/shani/NYU/fMRI_data/GLM_output/HRF_est/parameters/';

HRFfileName = [hrfDir 'SPM_HRF_params_' subj '_V1V2V3_restrict_21'];
load(HRFfileName) % this loads the variable bestfit
hrfParams.rdelay = bestfit.params(1);
hrfParams.udelay = bestfit.params(2);
hrfParams.udispersion= bestfit.params(3);
hrfParams.rdispersion = bestfit.params(4);

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





