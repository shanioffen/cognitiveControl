% 2008March31 - adapted runGLM to do the all subj 4 at a time leaving out 1

function runGLM_allSubj_by4s(expt)
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
model = 'corrTrial';
blur = 1; % set to 1 in order to spatially blur the data, to 0 to leave unblurred
split = 0; % going to be tricky to figure out how to choose the right HRF if split the data

% set the directory and start mrLR
clear view
clear global MLR
clear global mrDEFAULTS
baseDir = setBaseDir(expt,subj);
dataDir = [baseDir subj '_' expt];
cd(dataDir)
mrGlobals
view = newView('Volume');

% first set the stim files (which ones depend on model being tested)
view = stim2mrLR(view, model);

% if doing fixed effects, need to keep track of how many runs each sub did, 
% so you know which HRF to apply when running the GLM
% *** order of subjects is JG, DS, LM, RS, SO ********
% need to figure out what to do if split or run by 4s...
switch expt
  case 'detect'
    subjRunList = [10 17 10 10 20];
  case 'memory'
    subjRunList = [11 14 12 11 10];
end % switch expt
  

for iOut = 1:5 % go through, leaving out one subj at a time
  % set the junk frames and concatenate if haven't already
  if 1 
    view = junkAndCat(view, blur, split, iOut, subjRunList);
  end

  % keep only the runs for the subjects included
  params.subjRunList = subjRunList([1:(iOut-1) (iOut+1):5]);

  % set parameters to run GLM
  groupNames = viewGet(view,'groupNames');
  params.saveName = ['glm_spmHRF_nogui_' model];
  params.hrfModel = 'spmHRF_so';
  params.trSupersampling = 1;
  params.groupName = 'ConcatBlurByFours';
  params.whichModel = model;
  
  paramsInfo = {...
      {'groupName',groupNames,'Name of group from which to do eventRelated analysis'},...
      {'saveName',params.saveName,'File name to try to save as'},...
      {'hrfModel',params.hrfModel,'Name of the function that defines the hrf used in glm'},...
      {'trSupersampling', params.trSupersampling, 'minmax=[1 100]', 'TR supersampling factor (1=no supersampling) reulting design matrix will be downsampled afterwards'}};
  
  % load the hrf parameters for the subject (subfunction)
  params.hrfParams = loadHRFparams(subj,iOut);
  params.hrfParams.description = params.hrfModel;
  params.hrfParams.incDeriv = 0;
  params.hrfParams.paramInfo  = feval(params.hrfModel, 'params');
  
  view = viewSet(view,'groupName',params.groupName);
  groupNum = viewGet(view,'currentGroup');
  params.scanNum =viewGet(view,'nScans',groupNum); % default to last scan in group
  %  params.scanNum = 1:viewGet(view,'nScans',groupNum); % alternatively, take all scans in group
  
  useDefault = 1;
  params.scanParams = getEventRelatedParams(view,params,useDefault);
  
  for i = 1:length(params.scanNum)
    params.scanParams{params.scanNum(i)}.scanNum = params.scanNum(i);
    params.scanParams{params.scanNum(i)}.description = ['GLM analysis of ConcatBlur: ' num2str(i)];
    params.scanParams{params.scanNum(i)}.hdrlen = 30; % how far out to go, in seconds
  end
  
  % run the GLM
  [view d] = eventRelatedGlm_so(view,params);
  
end % for iOut = leaveOutList

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
function view = junkAndCat(view,blur, split, iOut, subjRunList)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% change junk frames and concatenate, leaving out one subj at a time

global MLR
for groupNum = viewGet(view,'groupNum','MotionComp');
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
params.newGroupName = 'ConcatBlurByFours'; 
params.description = 'Concatenation of [x...x] with hipass filtering'; 
params.filterType = 1; 
params.filterCutoff = 0.0100;
params.percentSignal = 1;
params.warp = 1; 
params.warpInterpMethod = 'linear';
params.blur = blur;

if split == 0 % concatenate all the data but leave out one subj at a time
  params.scanList = [1:sum(subjRunList(1:iOut-1)) sum(subjRunList(1:iOut))+1:sum(subjRunList)]+1;
  params.warpBaseScan = params.scanList(1); % choose first scan in list as base
  view = concatNblur(view,params);
else % split the data in half, every other run - not sure how to do this right yet
  for iSplit = 1:2
    params.scanList = 1+iSplit:2:(viewGet(view,'nScans',viewGet(view,'groupNum',params.groupName))-1);
    params.warpBaseScan = params.scanList(1); % choose first scan in list as base
    view = concatNblur(view,params);
  end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hrfParams = loadHRFparams(subj,iOut)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hrfDir = '/Users/shani/NYU/fMRI_data/GLM_output/HRF_est/parameters/';

if strcmp(subj,'allSubj')
  % have to deal separately with each subject
  subjList = {'JG','DS','LM','RS','SO'};% switched order because looking on JG brain so put first
  counter = 0;
  for iSub = [1:(iOut-1) (iOut+1):length(subjList)]
    counter = counter+1;
    sub = subjList{iSub};
    HRFfileName = [hrfDir 'SPM_HRF_params_' sub '_V1V2V3_restrict_21'];
    load(HRFfileName) % this loads the variable bestfit;
    hrfParams.rdelay(counter,1) = bestfit.params(1);
    hrfParams.udelay(counter,1) = bestfit.params(2);
    hrfParams.udispersion(counter,1) = bestfit.params(3);
    hrfParams.rdispersion(counter,1) = bestfit.params(4);
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

  % make sure we are running on a set with a stimfile
  stimfile = viewGet(view,'stimfile',params.scanNum(scanNum));
  
  if isempty(stimfile)
    mrMsgBox(sprintf('No associated stimfile with scan %i in group %s',params.scanNum(scanNum),params.groupName));
    scanParams = [];
    return
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % see if we have a stimfile from mgl, in which case we should
  % ask the user what the variable name is that they want ot use for the analysis
  if strfind(stimfile{1}.filetype,'mgl')

    % check to see what style this is, if the task variable does
    % not have a segmentTrace then it mus be an old style, in which
    % we used channels
    task = cellArray(stimfile{1}.task,2);
    if isfield(stimfile{1}.myscreen,'traces') && ~isfield(task{1}{1},'segmentTrace')
      % this is the old style, get the stimtrace number
      taskVarParams{end+1} = {'stimtrace',stimfile{1}.myscreen.stimtrace,'the trace number that contains the stimulus','incdec=[-1 1]',sprintf('minmax=[%i %i]',stimfile{1}.myscreen.stimtrace,size(stimfile{1}.myscreen.traces,1))};
    else
      % this is the new tyle, ask for a variable name
      [varnames varnamesStr] = getTaskVarnames(stimfile{1}.task);
      % if there is more than one task, then ask the user for that
      task = cellArray(stimfile{1}.task,2);
      if length(task)>1
	taskVarParams{end+1} = {'taskNum',num2cell(1:length(task)),'The task you want to use'};
      end
      % if there are multiple phases, then ask for that
      maxPhaseNum = 0;
      maxSegNum = 0;
      for tnum = 1:length(task)
	phaseNum{tnum} = num2cell(1:length(task{tnum}));
	maxPhaseNum = max(maxPhaseNum,length(task{tnum}));
	% if there are multiple _segments_, then ask for that
	for pnum = 1:length(task{tnum})
	  segNum{tnum}{pnum} = num2cell(1:length(task{tnum}{pnum}.segmin));
	  maxSegNum = max(maxSegNum,length(segNum{tnum}{pnum}));
	end
      end
      if maxPhaseNum > 1
	if length(task) == 1
	  taskVarParams{end+1} = {'phaseNum',phaseNum{1},'The phase of the task you want to use'};
	else
	  taskVarParams{end+1} = {'phaseNum',phaseNum,'The phase of the task you want to use','contingent=taskNum'};
	end
      end
      
      % if there is more than one segement in any of the phases, ask the user to specify
      % should add some error checking.
      if maxSegNum > 1
	  taskVarParams{end+1} = {'segmentNum',1,'The segment of the task you want to use','incdec=[-1 1]'};
      end
      
      % set up to get the variable name from the user
      taskVarParams{end+1} ={'varname',varnames{1},sprintf('Analysis variables: %s',varnamesStr)};
    end
  end %checking for mgl file
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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




