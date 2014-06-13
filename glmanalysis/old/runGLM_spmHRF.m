% 2007August19
% changed to have new HRF
% 2008March10 - expanded to use different models (Added the 'which')

function runGLM_spmHRF

exptList = {'memory','detect','memDetect','vertical'};
subjList = {'DS','JG','LM','RS','SO'};
total = length(subjList)*length(exptList);
whichModel = {'default','corrDelay','corrTrial'};
% (the difference between models:) 
% default uses 3 predictors: s1, delay, s2 (treats all trials the same)
% correctDelay uses 4 predictors: s1, delay_correct, delay_wrong, s2
% correctTrial uses 6 predictors: s1_correct, delay_correct, s2_correct, and s1,d,s2 _wrong
view = [];
MLR = [];
mrDEFAULTS = [];
for iExpt = 1:2 % 1:length(exptList)
  expt = exptList{iExpt};
  for iSubj = 1:5 % 1:length(subjList)
    subj = subjList{iSubj};
    for iModel = 2:3
      model = whichModel{iModel};
    
      % first set the stim files (which ones depend on model being tested)
      stim2mrLR(expt,subj,model);
    
      % note if the subject ran across 2 days (because then need to warp when concatenate)
      ran2 = 0;
      switch subj
        case 'DS', if strcmp(expt,'memory')||strcmp(expt,'detect'), ran2 = 1; end
        case 'RS', if strcmp(expt,'detect'), ran2 = 1; end
        case 'SO', if strcmp(expt,'detect'), ran2 = 1; end
      end
      
      % set the directory
      clear view MLR mrDEFAULTS
      baseDir = setBaseDir(expt);
      dataDir = [baseDir subj '_' expt];
      cd(dataDir)
      
      mrGlobals
      view = newView('Volume');
      
      % set the junk frames and concatenate if haven't already
      if(0)
      view = junkAndCat(view, ran2);
      end
      
      % set parameters to run GLM
      groupNames = viewGet(view,'groupNames');
      params.saveName = ['glm_spmHRF_nogui_' model];
      params.hrfModel = 'spmHRF_so';
      params.trSupersampling = 1;
      params.groupName = 'Concatenation';
      params.whichModel = model;
      
      paramsInfo = {...
          {'groupName',groupNames,'Name of group from which to do eventRelated analysis'},...
          {'saveName',params.saveName,'File name to try to save as'},...
          {'hrfModel',params.hrfModel,'Name of the function that defines the hrf used in glm'},...
          {'trSupersampling', params.trSupersampling, 'minmax=[1 100]', 'TR supersampling factor (1=no supersampling) reulting design matrix will be downsampled afterwards'}};
      
      % load the hrf parameters for the subject"
      hrfDir = '/Users/shani/NYU/fMRI_data/GLM_output/HRF_est/parameters/';
      HRFfileName = [hrfDir 'SPM_HRF_params_' subj '_V1V2V3_restrict_21'];
      load(HRFfileName) % this loads the variable bestfit
      
      params.hrfParams.description = params.hrfModel;
      params.hrfParams.rdelay = bestfit.params(1);
      params.hrfParams.udelay = bestfit.params(2);
      params.hrfParams.udispersion= bestfit.params(3);
      params.hrfParams.rdispersion = bestfit.params(4);
      params.hrfParams.incDeriv = 0;
      params.hrfParams.paramInfo  = feval(params.hrfModel, 'params');
      
      view = viewSet(view,'groupName',params.groupName);
      groupNum = viewGet(view,'currentGroup');
      params.scanNum =viewGet(view,'nScans',groupNum); % default to last scan in group 
      % params.scanNum = 1:viewGet(view,'nScans',groupNum); % alternatively, take all scans in group
      
      useDefault = 1;
      params.scanParams = getEventRelatedParams(view,params,useDefault);
      
      for i = 1:length(params.scanNum)
        params.scanParams{params.scanNum(i)}.scanNum = params.scanNum(i);
        params.scanParams{params.scanNum(i)}.description = ['GLM analysis of Concatenation: ' num2str(i)];
        params.scanParams{params.scanNum(i)}.hdrlen = 30; % how far out to go, in seconds
      end
      
      % run the GLM
      [view d] = eventRelatedGlm_so(view,params);
      
      
    end % iModel
  end % iSubj
end % iExpt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set directories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function baseDir = setBaseDir(expt)
switch expt
 case 'memory'
  baseDir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/';
 case 'detect'
  baseDir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/';
 case 'memDetect'
  baseDir = '/Users/shani/NYU/fMRI_data/MemDetect/';
 case 'vertical'
  baseDir = '/Users/shani/NYU/fMRI_data/Vertical/';
end % switch expt
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% change junk frames and concatenate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function view = junkAndCat(view,ran2)
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
params.newGroupName = 'Concatenation';
params.description = 'Concatenation of [x...x] with hipass filtering';
params.filterType = 1;
params.filterCutoff = 0.0100;
params.percentSignal = 1;
params.warp = 1; params.warpBaseScan = 2; 
params.warpInterpMethod = 'linear';
params.scanList = 2:(viewGet(view,'nScans',viewGet(view,'groupNum',params.groupName))-1);

view = concatTSeries(view,params);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to get the variable name that the user wants
% to do the event related analysis on, puts up a gui
% TAKEN FROM eventRelatedGlmGUI.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function scanParams = getEventRelatedParams(view,params,useDefault)

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




