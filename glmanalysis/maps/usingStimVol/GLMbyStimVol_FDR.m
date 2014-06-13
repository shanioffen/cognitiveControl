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

function [view d] = GLMbyStimVol_FDR(view,params)

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
params = mrParamsReconcile([],params);
if ieNotDefined('params'),return,end

% set the group
view = viewSet(view,'groupName',params.groupName);

% set the overlays
[r2 pval pvalTH correlation dpCorrect] = setOverlays(view,params); 

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
  ehdr = [];ehdrste = [];thisr2 = [];thisCorrelation = []; thisPval = [];

  for i = 1:ceil(numSlices/numSlicesAtATime)
    % load the scan
    d = loadScan(view,scanNum,[],[currentSlice min(numSlices,currentSlice+numSlicesAtATime-1)]);
    params.hrfParams.tmax = params.scanParams{scanNum}.hdrlen; %+d.tr/2;
    
    d.supersampling = params.trSupersampling;
    % use the duration of stimuli/events in the design matrix
    d.impulse = 0; 

    % get the stim volumes, if empty then abort
    d = getStimvol(d,params.scanParams{scanNum});
    if isempty(d.stimvol),mrWarnDlg('No stim volumes found');return,end
    
    % do any called-for preprocessing
    hrf = feval(params.hrfModel, d.tr/d.supersampling, params.hrfParams);
    d = eventRelatedPreProcess(d,params.scanParams{scanNum}.preprocess);
    % make a stimulation convolution matrix
    
    d = makeglm_so(d,hrf,params.subjRunList);
    % compute the betas
    d = getr2(d);
    % get contrast coeff so can get pvals
    d = getCorrCoef(d);    
    % update the current slice we are working on
    currentSlice = currentSlice+numSlicesAtATime;
    % cat with what has already been computed for other slices
    ehdr = cat(3,ehdr,d.ehdr);
    ehdrste = cat(3,ehdrste,d.ehdrste);
    thisr2 = cat(3,thisr2,d.r2);
    thisCorrelation = cat(3,thisCorrelation,d.correlation);
    thisPval = cat(3,thisPval,d.pval);
  end

  % now put all the data from all the slices into the structure
  d.ehdr = ehdr;
  d.ehdrste = ehdrste;
  d.r2 = thisr2;
  d.correlation = thisCorrelation;
  d.pval = thisPval;

  d.dim(3) = size(d.r2,3);

  % save the r2 overlay
  r2.data{scanNum} = d.r2;
  r2.params{scanNum} = params.scanParams{scanNum};
  
  % save the correlation overlay
  correlation.data{scanNum} = d.correlation;
  correlation.params{scanNum} = params.scanParams{scanNum};

  % set pval overlay
  pval.data{scanNum} = d.pval;
  pval.params{scanNum} = params.scanParams{scanNum};
  
  % get FDR threshold
  pvals = d.pval(:);
  [TH discard] = FDR(pvals,0.05); clear pvals 
  
  % calculate overlay using that TH 
  pvalTH.data{scanNum} = TH-pval.data{scanNum}; 
  pvalTH.params{scanNum} = params.scanParams{scanNum};  
  
  % also keep beta estimate for delay period
  dpCorrect.data{scanNum} = squeeze(d.ehdr(:,:,:,2));
  dpCorrect.params{scanNum} = params.scanParams{scanNum};
  
  
  
  
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
  erAnal.d{scanNum}.stimNames = d.stimNames;
  erAnal.d{scanNum}.scm = d.scm;
  erAnal.d{scanNum}.expname = d.expname;
  erAnal.d{scanNum}.fullpath = d.fullpath;
  erAnal.d{scanNum}.FDR_TH = TH;
  stimvol = d.stimvol;
  for i=1:length(stimvol)
      stimvol{i} = unique(ceil(stimvol{i}/d.supersampling));
  end
  erAnal.d{scanNum}.stimvol = stimvol;

end
toc

% set ranges
r2.range = [0 1];
correlation.range = [-1 1];
pval.range = [0 1];
pvalTH.range = [TH-1 TH];
dpCorrect.range = findRange(dpCorrect.data);


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
view = viewSet(view,'newoverlay',correlation);
view = viewSet(view,'newoverlay',pvalTH);
view = viewSet(view,'newoverlay',pval);
view = viewSet(view,'newoverlay',dpCorrect);

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
function [r2 pval pvalTH correlation dpCorrect] = setOverlays(view,params); 
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
r2.interrogator = 'plotGLMbyDM_FDR';
r2.mergeFunction = 'defaultMergeParams';

pval.name = 'pval';
pval.groupName = params.groupName;
pval.function = 'eventRelatedGlm_so';
pval.reconcileFunction = 'defaultReconcileParams';
pval.data = cell(1,viewGet(view,'nScans'));
pval.date = dateString;
pval.params = cell(1,viewGet(view,'nScans'));
pval.colormap = hot(312);
pval.colormap = pval.colormap(end-255:end,:);
pval.alpha = 1;
pval.colormapType = 'setRangeToMax';
pval.interrogator = 'plotGLMbyDM_FDR';
pval.mergeFunction = 'defaultMergeParams';

pvalTH.name = 'pvalTH';
pvalTH.groupName = params.groupName;
pvalTH.function = 'eventRelatedGlm_so';
pvalTH.reconcileFunction = 'defaultReconcileParams';
pvalTH.data = cell(1,viewGet(view,'nScans'));
pvalTH.date = dateString;
pvalTH.params = cell(1,viewGet(view,'nScans'));
pvalTH.colormap = hot(312);
pvalTH.colormap = pvalTH.colormap(end-255:end,:);
pvalTH.alpha = 1;
pvalTH.colormapType = 'setRangeToMax';
pvalTH.interrogator = 'plotGLMbyDM_FDR';
pvalTH.mergeFunction = 'defaultMergeParams';

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
dpCorrect.interrogator = 'plotGLMbyDM_FDR';
dpCorrect.mergeFunction = 'defaultMergeParams';
    
correlation.name = 'correlation';
correlation.groupName = params.groupName;
correlation.function = 'eventRelatedGlm_so';
correlation.reconcileFunction = 'defaultReconcileParams';
correlation.data = cell(1,viewGet(view,'nScans'));
correlation.date = dateString;
correlation.params = cell(1,viewGet(view,'nScans'));
correlation.colormap = hot(312);
correlation.colormap = correlation.colormap(end-255:end,:);
correlation.alpha = 1;
correlation.colormapType = 'setRangeToMax';
correlation.interrogator = 'plotGLMbyDM_FDR';
correlation.mergeFunction = 'defaultMergeParams';

