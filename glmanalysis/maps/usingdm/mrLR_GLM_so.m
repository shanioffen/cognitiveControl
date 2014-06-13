% eventRelatedGlm.m
%
%      usage: view = eventRelatedGlm_so(view,params)
%         by: shani, modeled on Glm code by farshad moradi, modeled on eventRelated code by Justin
%       date: 2007August3
%    purpose: farshad's GLM is same as eventRelated, but uses canonical hrf instead of deconvolution
%             this code is the same as GLM, but uses measured HRF instead of canonical, and outputs
%             overlays of the beta values in addition to the rSquared overlay
%             TODO: need to also make it possible to input ROIs and do the GLM on the ROIs as well as voxels
% 2007Aug19 - saved a backup (eventRelatedGlm_so.BKP) and started editing to use SPM HRF
% 2008Mar10 - made it possible to run different models
% 2008Mar19 - made it possible to do fixed effects
% 2008Apr07 - using my own DM instead of making it Farshad's way...

function [view d] = mrLR_GLM_so(view,params)

d = [];

% check arguments
if ~any(nargin == [1 2])
  help eventRelated
  return
end

mrGlobals;

% Default parameters if not input by user:
if ieNotDefined('params')
  % put up the gui
  params = eventRelatedGlmGUI;
  params.whichModel = 'default'; % default to treat all trials the same;
  params.subjRunList = []; % default is just one subject
end

% Reconcile params and abort if empty
% params = defaultReconcileParams([],params); ******** NEED TO DEBUG AND PUT THIS BACK IN
if ieNotDefined('params'),return,end

% set the group
view = viewSet(view,'groupName',params.groupName);

% set the overlays
% now depends on the model being run
model = params.whichModel;
switch model
  case 'default'
    [r2 dpContrast s1 dp s2] = setOverlays(view,params); 
  case {'correctDelay','correctDelayS2plus','correctDelayS2minus'}
    [r2 dpContrast s1 dpCorrect dpWrong s2] = setOverlays(view,params); 
  case {'correctTrial','correctTrialS2plus','correctTrialS2minus'}
    [r2 dpContrast s1Correct dpCorrect s2Correct s1Wrong dpWrong s2Wrong] = setOverlays(view,params); 
end

tic
set(viewGet(view,'figNum'),'Pointer','watch');drawnow;
for scanNum = params.scanNum
  % decide how many slices to do at a time, this is done
  % simply to save memory -- currently our system is limited
  % to 2G of memory and for large concatenations, you need
  % to break up the analysis into smaller portions of the data
  numSlices = viewGet(view,'nSlices',scanNum);
  numVolumes = viewGet(view,'nFrames',scanNum);
  dims = viewGet(view,'dims',scanNum);
  % choose how many slices based on trying to keep a certain
  % amount of data in the memory
  numSlicesAtATime = getNumSlicesAtATime(numVolumes,dims);
  currentSlice = 1;
  ehdr = [];ehdrste = [];thisr2 = [];

  for i = 1:ceil(numSlices/numSlicesAtATime)
    % load the scan
    d = loadScan_noStimFile(view,scanNum,[],[currentSlice min(numSlices,currentSlice+numSlicesAtATime-1)]);
    params.hrfParams.tmax = params.scanParams{scanNum}.hdrlen; %+d.tr/2;
    
    d.supersampling = params.trSupersampling;
    % use the duration of stimuli/events in the design matrix
    d.impulse = 0; 

    % set DM from params
    d.DM = params.DM;
    
    % do any called-for preprocessing
    hrf = feval(params.hrfModel, d.tr/d.supersampling, params.hrfParams);
    d = eventRelatedPreProcess(d,params.scanParams{scanNum}.preprocess);
    
    % compute the glm
    disp('(mrLR_GLM_so) calling makeglmByDM_so')     
    d = makeglmByDM_so(d,hrf,params.subjRunList,scanNum);
    % compute the estimated hemodynamic responses
    d = getr2(d);
    % update the current slice we are working on
    currentSlice = currentSlice+numSlicesAtATime;
    % cat with what has already been computed for other slices
    ehdr = cat(3,ehdr,d.ehdr);
    ehdrste = cat(3,ehdrste,d.ehdrste);
    thisr2 = cat(3,thisr2,d.r2);
  end

  % now put all the data from all the slices into the structure
  d.ehdr = ehdr;
  d.ehdrste = ehdrste;
  d.r2 = thisr2;

  d.dim(3) = size(d.r2,3);

  % save the r2 overlay
  r2.data{scanNum} = d.r2;
  r2.params{scanNum} = params.scanParams{scanNum};
  
  % save the beta overlays; depends on model
  switch model
    case 'default'
      s1.data{scanNum} = squeeze(d.ehdr(:,:,:,1));
      s1.params{scanNum} = params.scanParams{scanNum};
      
      dp.data{scanNum} = squeeze(d.ehdr(:,:,:,2));
      dp.params{scanNum} = params.scanParams{scanNum};
      
      s2.data{scanNum} = squeeze(d.ehdr(:,:,:,3));
      s2.params{scanNum} = params.scanParams{scanNum};
      
    case {'correctDelay','correctDelayS2plus','correctDelayS2minus'}
      s1.data{scanNum} = squeeze(d.ehdr(:,:,:,1));
      s1.params{scanNum} = params.scanParams{scanNum};
      
      dpCorrect.data{scanNum} = squeeze(d.ehdr(:,:,:,2));
      dpCorrect.params{scanNum} = params.scanParams{scanNum};
      
      dpWrong.data{scanNum} = squeeze(d.ehdr(:,:,:,3));
      dpWrong.params{scanNum} = params.scanParams{scanNum};
      
      s2.data{scanNum} = squeeze(d.ehdr(:,:,:,4));
      s2.params{scanNum} = params.scanParams{scanNum};
      
      dpContrast.data{scanNum} = dpCorrect.data{scanNum} - dpWrong.data{scanNum};
      dpContrast.params{scanNum} = params.scanParams{scanNum};
      
    case {'correctTrial','correctTrialS2plus','correctTrialS2minus'}
      s1Correct.data{scanNum} = squeeze(d.ehdr(:,:,:,1));
      s1Correct.params{scanNum} = params.scanParams{scanNum};
      
      dpCorrect.data{scanNum} = squeeze(d.ehdr(:,:,:,2));
      dpCorrect.params{scanNum} = params.scanParams{scanNum};
      
      s2Correct.data{scanNum} = squeeze(d.ehdr(:,:,:,3));
      s2Correct.params{scanNum} = params.scanParams{scanNum};
            
      s1Wrong.data{scanNum} = squeeze(d.ehdr(:,:,:,4));
      s1Wrong.params{scanNum} = params.scanParams{scanNum};
      
      dpWrong.data{scanNum} = squeeze(d.ehdr(:,:,:,5));
      dpWrong.params{scanNum} = params.scanParams{scanNum};
      
      s2Wrong.data{scanNum} = squeeze(d.ehdr(:,:,:,6));
      s2Wrong.params{scanNum} = params.scanParams{scanNum};
      
      dpContrast.data{scanNum} = dpCorrect.data{scanNum} - dpWrong.data{scanNum};
      dpContrast.params{scanNum} = params.scanParams{scanNum};
      
  end %switch model
  
  % save other eventRelated parameters
  erAnal.d{scanNum}.hrf = d.simulatedhrf;
  erAnal.d{scanNum}.actualhrf = hrf;
  erAnal.d{scanNum}.trsupersampling = d.supersampling;
  erAnal.d{scanNum}.ver = d.ver;
  erAnal.d{scanNum}.filename = d.filename;
  erAnal.d{scanNum}.filepath = d.filepath;
  erAnal.d{scanNum}.dim = d.dim;
  erAnal.d{scanNum}.ehdr = d.ehdr;
  erAnal.d{scanNum}.ehdrste = d.ehdrste;
  erAnal.d{scanNum}.nhdr = d.nhdr;
  erAnal.d{scanNum}.hdrlen = d.hdrlen;
  erAnal.d{scanNum}.tr = d.tr;
  erAnal.d{scanNum}.stimNames = 1:size(d.DM,2);
  erAnal.d{scanNum}.scm = d.scm;
  erAnal.d{scanNum}.expname = d.expname;
  erAnal.d{scanNum}.fullpath = d.fullpath;
  erAnal.d{scanNum}.DM = d.DM;

end
toc

% set ranges
r2.range = [0 1];

% others depend on model
switch model
  case 'default'
    s1.range = findRange(s1.data);
    dp.range = findRange(dp.data);
    s2.range = findRange(s2.data);
  case {'correctDelay','correctDelayS2plus','correctDelayS2minus'}
    s1.range = findRange(s1.data);
    dpCorrect.range = findRange(dpCorrect.data);
    dpWrong.range = findRange(dpWrong.data);
    s2.range = findRange(s2.data);
    dpContrast.range = findRange(dpContrast.data);
  case {'correctTrial','correctTrialS2plus','correctTrialS2minus'}
    s1Correct.range = findRange(s1Correct.data);
    dpCorrect.range = findRange(dpCorrect.data);
    s2Correct.range = findRange(s2Correct.data);
    s1Wrong.range = findRange(s1Wrong.data);
    dpWrong.range = findRange(dpWrong.data);
    s2Wrong.range = findRange(dpWrong.data);
    dpContrast.range = findRange(dpContrast.data);    
end


% install analysis
erAnal.name = params.saveName;
erAnal.type = 'glmAnal';
erAnal.groupName = params.groupName;
erAnal.function = 'eventRelatedGlm';
erAnal.reconcileFunction = 'defaultReconcileParams';
erAnal.mergeFunction = 'defaultMergeParams';
erAnal.guiFunction = 'eventRelatedGlmGUI';
erAnal.params = params;
erAnal.date = r2.date;

view = viewSet(view,'newAnalysis',erAnal);
view = viewSet(view,'newoverlay',r2);
% set overlays depending on model
switch model
  case 'default'
    view = viewSet(view,'newoverlay',s1);
    view = viewSet(view,'newoverlay',dp);
    view = viewSet(view,'newoverlay',s2);
  case {'correctDelay','correctDelayS2plus','correctDelayS2minus'}
    view = viewSet(view,'newoverlay',s1);
    view = viewSet(view,'newoverlay',dpCorrect);
    view = viewSet(view,'newoverlay',dpWrong);
    view = viewSet(view,'newoverlay',s2);
    view = viewSet(view,'newoverlay',dpContrast);    
  case {'correctTrial','correctTrialS2plus','correctTrialS2minus'}
    view = viewSet(view,'newoverlay',s1Correct);
    view = viewSet(view,'newoverlay',dpCorrect);
    view = viewSet(view,'newoverlay',s2Correct);
    view = viewSet(view,'newoverlay',s1Wrong);
    view = viewSet(view,'newoverlay',dpWrong);
    view = viewSet(view,'newoverlay',s2Wrong);
    view = viewSet(view,'newoverlay',dpContrast);        
end    

% Save it
saveAnalysis(view,erAnal.name);

set(viewGet(view,'figNum'),'Pointer','arrow');drawnow

% for output
if nargout > 1
  for i = 1:length(d)
    erAnal.d{i}.r2 = r2.data{i};
  end
  % make d strucutre
  if length(erAnal.d) == 1
    d = erAnal.d{1};
  else
    d = erAnal.d;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         SUBFUNCTIONS                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function range = findRange(data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ampMin = realmax;
ampMax = 0;
nScans = length(data);
for scan=1:nScans
  if ~isempty(data{scan})
    ampMin = min([ampMin min(data{scan}(:))]);
    ampMax = max([ampMax max(data{scan}(:))]);
  end
end
if (ampMin <= ampMax)
  range = [ampMin ampMax];
else
  % if amp data is empty, need to make sure min < max
  range = [0 1];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [r2 dpContrast varargout] = setOverlays(view,params) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create the parameters for the r2 overlay
dateString = datestr(now);
r2.name = 'r2';
r2.groupName = params.groupName;
r2.function = 'eventRelatedGlm_so';
r2.reconcileFunction = 'defaultReconcileParams';
r2.data = cell(1,viewGet(view,'nScans'));
r2.date = dateString;
r2.params = cell(1,viewGet(view,'nScans'));
r2.colormap = hot(312); % colormap is made with a little bit less on the dark end
r2.colormap = r2.colormap(end-255:end,:);
r2.alpha = 1;
r2.colormapType = 'setRangeToMax';
r2.interrogator = 'plotGLMbyDM_so';
r2.mergeFunction = 'defaultMergeParams';

dpContrast.name = 'dpContrast';
dpContrast.groupName = params.groupName;
dpContrast.function = 'eventRelatedGlm_so';
dpContrast.reconcileFunction = 'defaultReconcileParams';
dpContrast.data = cell(1,viewGet(view,'nScans'));
dpContrast.date = dateString;
dpContrast.params = cell(1,viewGet(view,'nScans'));
dpContrast.colormap = hot(312);
dpContrast.colormap = dpContrast.colormap(end-255:end,:);
dpContrast.alpha = 1;
dpContrast.colormapType = 'setRangeToMax';
dpContrast.interrogator = 'plotGLMbyDM_so';
dpContrast.mergeFunction = 'defaultMergeParams';
    
% create the parameters for the beta overlays
% how many and which depend on the model being used
model = params.whichModel;
switch model
  case 'default'
    s1.name = 's1';
    s1.groupName = params.groupName;
    s1.function = 'eventRelatedGlm_so';
    s1.reconcileFunction = 'defaultReconcileParams';
    s1.data = cell(1,viewGet(view,'nScans'));
    s1.date = dateString;
    s1.params = cell(1,viewGet(view,'nScans'));
    s1.colormap = hot(312);
    s1.colormap = s1.colormap(end-255:end,:);
    s1.alpha = 1;
    s1.colormapType = 'setRangeToMax';
    s1.interrogator = 'plotGLMbyDM_so';
    s1.mergeFunction = 'defaultMergeParams';
    varargout{1} = s1;
    
    dp.name = 'dp';
    dp.groupName = params.groupName;
    dp.function = 'eventRelatedGlm_so';
    dp.reconcileFunction = 'defaultReconcileParams';
    dp.data = cell(1,viewGet(view,'nScans'));
    dp.date = dateString;
    dp.params = cell(1,viewGet(view,'nScans'));
    dp.colormap = hot(312);
    dp.colormap = dp.colormap(end-255:end,:);
    dp.alpha = 1;
    dp.colormapType = 'setRangeToMax';
    dp.interrogator = 'plotGLMbyDM_so';
    dp.mergeFunction = 'defaultMergeParams';
    varargout{2} = dp;
    
    s2.name = 's2';
    s2.groupName = params.groupName;
    s2.function = 'eventRelatedGlm_so';
    s2.reconcileFunction = 'defaultReconcileParams';
    s2.data = cell(1,viewGet(view,'nScans'));
    s2.date = dateString;
    s2.params = cell(1,viewGet(view,'nScans'));
    s2.colormap = hot(312);
    s2.colormap = s2.colormap(end-255:end,:);
    s2.alpha = 1;
    s2.colormapType = 'setRangeToMax';
    s2.interrogator = 'plotGLMbyDM_so';
    s2.mergeFunction = 'defaultMergeParams';
    varargout{3} = s2;
    
  case {'correctDelay','correctDelayS2plus','correctDelayS2minus'}
    s1.name = 's1';
    s1.groupName = params.groupName;
    s1.function = 'eventRelatedGlm_so';
    s1.reconcileFunction = 'defaultReconcileParams';
    s1.data = cell(1,viewGet(view,'nScans'));
    s1.date = dateString;
    s1.params = cell(1,viewGet(view,'nScans'));
    s1.colormap = hot(312);
    s1.colormap = s1.colormap(end-255:end,:);
    s1.alpha = 1;
    s1.colormapType = 'setRangeToMax';
    s1.interrogator = 'plotGLMbyDM_so';
    s1.mergeFunction = 'defaultMergeParams';
    varargout{1} = s1;
    
    dpCorrect.name = 'dpCorrect';
    dpCorrect.groupName = params.groupName;
    dpCorrect.function = 'eventRelatedGlm_so';
    dpCorrect.reconcileFunction = 'defaultReconcileParams';
    dpCorrect.data = cell(1,viewGet(view,'nScans'));
    dpCorrect.date = dateString;
    dpCorrect.params = cell(1,viewGet(view,'nScans'));
    dpCorrect.colormap = hot(312);
    dpCorrect.colormap = dpCorrect.colormap(end-255:end,:);
    dpCorrect.alpha = 1;
    dpCorrect.colormapType = 'setRangeToMax';
    dpCorrect.interrogator = 'plotGLMbyDM_so';
    dpCorrect.mergeFunction = 'defaultMergeParams';
    varargout{2} = dpCorrect;
    
    dpWrong.name = 'dpWrong';
    dpWrong.groupName = params.groupName;
    dpWrong.function = 'eventRelatedGlm_so';
    dpWrong.reconcileFunction = 'defaultReconcileParams';
    dpWrong.data = cell(1,viewGet(view,'nScans'));
    dpWrong.date = dateString;
    dpWrong.params = cell(1,viewGet(view,'nScans'));
    dpWrong.colormap = hot(312);
    dpWrong.colormap = dpWrong.colormap(end-255:end,:);
    dpWrong.alpha = 1;
    dpWrong.colormapType = 'setRangeToMax';
    dpWrong.interrogator = 'plotGLMbyDM_so';
    dpWrong.mergeFunction = 'defaultMergeParams';
    varargout{3} = dpWrong;
    
    s2.name = 's2';
    s2.groupName = params.groupName;
    s2.function = 'eventRelatedGlm_so';
    s2.reconcileFunction = 'defaultReconcileParams';
    s2.data = cell(1,viewGet(view,'nScans'));
    s2.date = dateString;
    s2.params = cell(1,viewGet(view,'nScans'));
    s2.colormap = hot(312);
    s2.colormap = s2.colormap(end-255:end,:);
    s2.alpha = 1;
    s2.colormapType = 'setRangeToMax';
    s2.interrogator = 'plotGLMbyDM_so';
    s2.mergeFunction = 'defaultMergeParams';
    varargout{4} = s2;
    
  case {'correctTrial','correctTrialS2plus','correctTrialS2minus'}
    s1Correct.name = 's1Correct';
    s1Correct.groupName = params.groupName;
    s1Correct.function = 'eventRelatedGlm_so';
    s1Correct.reconcileFunction = 'defaultReconcileParams';
    s1Correct.data = cell(1,viewGet(view,'nScans'));
    s1Correct.date = dateString;
    s1Correct.params = cell(1,viewGet(view,'nScans'));
    s1Correct.colormap = hot(312);
    s1Correct.colormap = s1Correct.colormap(end-255:end,:);
    s1Correct.alpha = 1;
    s1Correct.colormapType = 'setRangeToMax';
    s1Correct.interrogator = 'plotGLMbyDM_so';
    s1Correct.mergeFunction = 'defaultMergeParams';
    varargout{1} = s1Correct;
    
    dpCorrect.name = 'dpCorrect';
    dpCorrect.groupName = params.groupName;
    dpCorrect.function = 'eventRelatedGlm_so';
    dpCorrect.reconcileFunction = 'defaultReconcileParams';
    dpCorrect.data = cell(1,viewGet(view,'nScans'));
    dpCorrect.date = dateString;
    dpCorrect.params = cell(1,viewGet(view,'nScans'));
    dpCorrect.colormap = hot(312);
    dpCorrect.colormap = dpCorrect.colormap(end-255:end,:);
    dpCorrect.alpha = 1;
    dpCorrect.colormapType = 'setRangeToMax';
    dpCorrect.interrogator = 'plotGLMbyDM_so';
    dpCorrect.mergeFunction = 'defaultMergeParams';
    varargout{2} = dpCorrect;
    
    s2Correct.name = 's2Correct';
    s2Correct.groupName = params.groupName;
    s2Correct.function = 'eventRelatedGlm_so';
    s2Correct.reconcileFunction = 'defaultReconcileParams';
    s2Correct.data = cell(1,viewGet(view,'nScans'));
    s2Correct.date = dateString;
    s2Correct.params = cell(1,viewGet(view,'nScans'));
    s2Correct.colormap = hot(312);
    s2Correct.colormap = s2Correct.colormap(end-255:end,:);
    s2Correct.alpha = 1;
    s2Correct.colormapType = 'setRangeToMax';
    s2Correct.interrogator = 'plotGLMbyDM_so';
    s2Correct.mergeFunction = 'defaultMergeParams';
    varargout{3} = s2Correct;
    
    s1Wrong.name = 's1Wrong';
    s1Wrong.groupName = params.groupName;
    s1Wrong.function = 'eventRelatedGlm_so';
    s1Wrong.reconcileFunction = 'defaultReconcileParams';
    s1Wrong.data = cell(1,viewGet(view,'nScans'));
    s1Wrong.date = dateString;
    s1Wrong.params = cell(1,viewGet(view,'nScans'));
    s1Wrong.colormap = hot(312);
    s1Wrong.colormap = s1Wrong.colormap(end-255:end,:);
    s1Wrong.alpha = 1;
    s1Wrong.colormapType = 'setRangeToMax';
    s1Wrong.interrogator = 'plotGLMbyDM_so';
    s1Wrong.mergeFunction = 'defaultMergeParams';
    varargout{4} = s1Wrong;
    
    dpWrong.name = 'dpWrong';
    dpWrong.groupName = params.groupName;
    dpWrong.function = 'eventRelatedGlm_so';
    dpWrong.reconcileFunction = 'defaultReconcileParams';
    dpWrong.data = cell(1,viewGet(view,'nScans'));
    dpWrong.date = dateString;
    dpWrong.params = cell(1,viewGet(view,'nScans'));
    dpWrong.colormap = hot(312);
    dpWrong.colormap = dpWrong.colormap(end-255:end,:);
    dpWrong.alpha = 1;
    dpWrong.colormapType = 'setRangeToMax';
    dpWrong.interrogator = 'plotGLMbyDM_so';
    dpWrong.mergeFunction = 'defaultMergeParams';
    varargout{5} = dpWrong;
    
    s2Wrong.name = 's2Wrong';
    s2Wrong.groupName = params.groupName;
    s2Wrong.function = 'eventRelatedGlm_so';
    s2Wrong.reconcileFunction = 'defaultReconcileParams';
    s2Wrong.data = cell(1,viewGet(view,'nScans'));
    s2Wrong.date = dateString;
    s2Wrong.params = cell(1,viewGet(view,'nScans'));
    s2Wrong.colormap = hot(312);
    s2Wrong.colormap = s2Wrong.colormap(end-255:end,:);
    s2Wrong.alpha = 1;
    s2Wrong.colormapType = 'setRangeToMax';
    s2Wrong.interrogator = 'plotGLMbyDM_so';
    s2Wrong.mergeFunction = 'defaultMergeParams';
    varargout{6} = s2Wrong;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function hrf = spmHRF_so(TR,params) %

% **** now this is done in a separate function, not a subfunction, so doesn't get confusing since call it from many different functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

