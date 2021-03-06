function view = corrByVox(view,seed,params)
%
% view = corrByVox(view,seed,[params])
% 
% Loops throughs scans, loads corresponding tSeries, and computes,
% for each voxel, the correlation of its time series to the seed
%
% The seed can be either 
%    [1] a voxel (in volume space, [x y z]), ** not yet implemented because may not be useful **
% or [2] an ROI name (e.g. 'V1', must be a string of characters),
% or [3] a timecourse vector, in which case you must give a cell array with 
% seed{1}='timecourse name' (so it can be identified) and seed{2} = a timecourse vector 
% (of the correct length). ** It is assumed that this TC already has been processed in any
% way that is desired (e.g. FilterF detrending, divide by mean, etc) so nothing more will
% be done to it, unless you set params.detrendSeed = 1, then the same filter that was applied
% to the data (if any) when it was concatenated will also be applied to the seed time-course**
%
% params: Optional initial parameters. Default: current group, all scans, no filtering, subj JG
% params.groupName: group of scans that will be analyzed.
%                   Default: current group of view.
% params.scanList: Which scan to compute (integer index into scans)
%                   Default: last scan (n)
%
% params.detrendSeed: Now just inherit whatever detrending was done on the scans being analyzed (Oct 2007)
%                     for an ROI or a voxel, don't need to do anything, because pulling the seed from the same scans as the data
%                     for a time-course seed, can set this == 1, and then use the same filter as was used when
%                     the data were extracted by concatenation
%
% params.subj:       Need to input the subj so can load the proper HRF for creating the DM for plotting
%
% shani, August 2007, Oct 2007

%**********************************************
% I'm not making this to run from within the GUI, so 
% there needs to be some sort of wrapper that calls MLR etc
% may make this separate code, but meanwhile, it'll be something like this
mrGlobals

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get analysis parameters

% if no subject is defined, default to JG for now
if ieNotDefined('params.subj'), params.subj = 'JG'; end

% define the hrf parameters for the subject
params.hrfModel = 'spmHRF_so';
params.trSupersampling = 1;
hrfDir = '~/NYU/fMRI_data/GLM_output/HRF_est/parameters/';
HRFfileName = [hrfDir 'SPM_HRF_params_' params.subj '_V1V2V3_restrict_21'];
load(HRFfileName) % this loads the variable bestfit
params.hrfParams.description = params.hrfModel;
params.hrfParams.rdelay = bestfit.params(1);
params.hrfParams.udelay = bestfit.params(2);
params.hrfParams.udispersion= bestfit.params(3);
params.hrfParams.rdispersion = bestfit.params(4);
params.hrfParams.incDeriv = 0;
params.hrfParams.paramInfo  = feval(params.hrfModel, 'params');
clear bestfit

% default to current group:
if ieNotDefined('params.groupName'), params.groupName = viewGet(view,'groupName'); end

% default to last scan in group:
if ieNotDefined('params.scanList') 
  groupNum = viewGet(view,'groupNum',params.groupName); 
  params.scanList = viewGet(view,'nScans',groupNum); 
end

% Change group
groupName = params.groupName;
curGroup = viewGet(view,'currentGroup');
groupNum = viewGet(view,'groupNum',groupName);
if (groupNum ~= curGroup)
  mrWarnDlg(['Changing view to group: ',groupName]);
  view = viewSet(view,'currentGroup',groupNum);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check what kind of seed we have and set the seed time course accordingly:

if ischar(seed)
  % then the user input an ROI, and we need to extract the mean time course for that ROI:
  seedType = ['ROI:' seed];			
  seedTC = getROITimeCourse(view,seed,params); 
elseif length(seed) == 3
  % then the user has input a voxel, and we need to extract the time course for that voxel
  seedType = ['Voxel:' num2str(seed(1)) '-' num2str(seed(2)) '-' num2str(seed(3))];
  seedTC = loadVoxTimeCourse(view,seed,params); 
  display('I havent implemented a voxel seed yet; seed must be an ROI or a timecourse')
  return
elseif iscell(seed)
  seedType = seed{1};
  seedTC = makeTCcell(view, seed{2}, params);
end

% add the seed to the parameters
params.seedTC = seedTC; % this is a cell array of the same length as scanList
params.seedType = seedType;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create the analysis structure
dateString = datestr(now);
corByVox.name = ['CorrToSeed_' seedType];  % This can be reset by editAnalysisGUI
corByVox.type = 'CorrelationToSeed';
corByVox.function = 'corrByVox';
corByVox.groupName = params.groupName;
corByVox.guiFunction = '';
corByVox.params = params;
corByVox.date = dateString;
[tf corByVox] = isanalysis(corByVox); % in case need to add fields

% need to create a data structure with the design matrix so can do 
% trial-triggered-average plots with the interrogator
corByVox.d = mkD(view, params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create the overlay
corVal.name = 'corVal';
corVal.function = 'corrbyVox';
corVal.date = dateString;
corVal.params = params;
corVal.colormap = jet(256);
corVal.interrogator = 'corPlot_so';
corVal.groupName = params.groupName;
corVal.data = cell(1,max(params.scanList)); % initialize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For each voxel in each scan in the group, calculate
% its correlation to the seed
disp('Computing correlation to seed ...');
waitHandle = mrWaitBar(0,'Computing correlation to seed');
warning('off','MATLAB:divideByZero');

for iScan = 1:length(params.scanList) 
  scanNum = params.scanList(iScan)
  params.concatInfo = viewGet(view,'concatInfo',scanNum);
  
  % sliceDims: [ydim xdim] for single slice
  % volDims; [ydim xdim nslices] for single scan
  sliceDims = viewGet(view,'sliceDims',scanNum);
  volDims = viewGet(view,'dims',scanNum);
  
  % Initialize data with NaNs
  corVal.data{scanNum} = NaN*ones(volDims);
  
  nslices = viewGet(view,'nslices',scanNum);
  waitTotal = nslices*length(params.scanList); waitCount = 0;
  for sliceNum = 1:nslices
    waitCount = waitCount + 1;
    temp = computeCorByVox(view,scanNum,iScan,sliceNum,params);

    corVal.data{scanNum}(:,:,sliceNum) = reshape(temp,sliceDims);
    clear temp
    mrWaitBar(waitCount/waitTotal,waitHandle);
    
  end % going through slices
end % going through scans

% set overlay range
corVal.range = findRange(corVal.data);


% Install analysis in the view
view = viewSet(view,'newanalysis',corByVox);
view = viewSet(view,'newoverlay',corVal);

% Save it
saveAnalysis(view,corByVox.name);

close(waitHandle)
return;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               SUBFUNCTIONS                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function seedTC = getROITimeCourse(view,seed,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display(['getting average TC for ROI ' seed])
roiName = seed;  % [seed '_restrict'];
view = loadROI(view,roiName);
roiNum = viewGet(view,'roiNum',roiName); 
groupNum = viewGet(view,'groupNum',params.groupName);

temp = meanTSeries(view,groupNum,roiNum,params.scanList,'detrend','None');
seedTC = temp(1,:); clear temp % just to get it to the right shape
% note - don't need to detrend, just inherit the detrending from the concatenated data, 
% because just pulling the ROI
% By default, meanTSeries will divide by mean and subtract the mean

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function seedTC = loadVoxTimeCourse(view,seed,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

seedTC = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function seedTC = makeTCcell(view, TC, params)

if(~isfield(params,'detrendSeed'))
  params.detrendSeed = 0;
end

% turn a TC into a cell of TCs for each scan
% (expectation is that the user has already done any
% desired processing on the input TC, so no more is done)
seedTC = cell(1,length(params.scanList));
for ind = 1:length(params.scanList)
  
  % make sure seed is the right length:
  if length(TC) ~= viewGet(view,'nFrames',ind)
    error('Seed is not the right length')
  end

  % can choose to filter this seed in the same way as the data was detrended
  if(params.detrendSeed==1)
    concatInfo = viewGet(view,'concatInfo',params.scanList(ind));
    hpf = concatInfo.hipassfilter; 
    numRuns = length(hpf); numTP = length(hpf{1});
    for iRun = 1:numRuns
      filter(((iRun-1)*numTP)+1:(iRun*numTP),1) = [hpf{iRun}'];
    end
    seedTC{ind} = real(ifft(fft(TC) .* filter));
  else
    seedTC{ind} = TC;
  end
  
end


 
% NOTE: the seed time course is arranged by the index of scans,
% not by the scan number (eg, if there are 12 scans to choose from,
% and the scanList chooses scans 3,8 and 9, then the seedTC will 
% have three cells, the first for scan 3, the second for scan 8 and 
% the third for scan 9). This is in contrast to the analysis overlay,
% corVal.data, which will have blank fields for the scans which are excluded.
% In the example here, it will have 9 cells, and only the 3rd, 8th and 9th
% will actually have data values.

% This is because the function meanTseries returns the tSeries by
% index, not by scan number, so the seedTC must do it that way; however
% the analysis overlays must keep track of which scan the overlay is for,
% so it must have the extra cells.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d_out = mkD(view, params)

for iScan = 1:length(params.scanList) 
  scanNum = params.scanList(iScan);

  % get initial scan info for the data structure, but don't actually load the data
  d = loadScan(view,scanNum,[],0);
  params.hrfParams.tmax = 30; %+d.tr/2;
    
  d.supersampling = params.trSupersampling;
  % use the duration of stimuli/events in the design matrix
  d.impulse = 0; 

  % get the stim volumes, if empty then abort
  d = getStimvol(d,scanNum);
  if isempty(d.stimvol),mrWarnDlg('No stim volumes found');return,end
    
  % do any called-for preprocessing
  hrf = feval(params.hrfModel, d.tr/d.supersampling, params.hrfParams);
  
  % make the design matrix
  d = makeglm(d,hrf);
  
  % save other eventRelated parameters
  d_out{scanNum}.hrf = d.simulatedhrf;
  d_out{scanNum}.actualhrf = hrf;
  d_out{scanNum}.trsupersampling = d.supersampling;
  d_out{scanNum}.ver = d.ver;
  d_out{scanNum}.filename = d.filename;
  d_out{scanNum}.filepath = d.filepath;
  d_out{scanNum}.dim = d.dim;
  d_out{scanNum}.tr = d.tr;
  d_out{scanNum}.stimNames = d.stimNames;
  d_out{scanNum}.scm = d.scm;
  d_out{scanNum}.expname = d.expname;
  d_out{scanNum}.fullpath = d.fullpath;

  stimvol = d.stimvol;
  for i=1:length(stimvol)
      stimvol{i} = unique(ceil(stimvol{i}/d.supersampling));
  end
  d_out{scanNum}.stimvol = stimvol;

end % going through scans


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function corVal = computeCorByVox(view, scan, iScan, slice, params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computes correlation to the seed for the given scan and slice.

% Analysis parameters for this scan
junkframes = viewGet(view,'junkframes',scan);
nframes = viewGet(view,'nframes',scan);

% Load tSeries
tSeries = loadTSeries(view, scan, slice);

% Reshape the tSeries
tSeries = reshapeTSeries(tSeries);

% Remove junkFrames
tSeries = tSeries(junkframes+1:junkframes+nframes,:);

% Remove dc, convert to percent, detrend, and spatial normalization
warnState = warning('query','MATLAB:divideByZero');
warning('off','MATLAB:divideByZero');

% the default when extracting the time series is to do spatial
% normalization by dividing by the mean; we leave that as is
% but we don't detrend because it's already been detrended when it was concatenated
ptSeries = percentTSeries(tSeries,...
    'detrend', 'None',...
    'subtractMean', 'Yes',...
    'temporalNormalization', 'No');
warning(warnState.state,warnState.identifier);
clear tSeries

% loop through to get correlation of each voxel to the seed 
% (faster than computing the correlation for the whole matrix)
for iVox = 1:size(ptSeries,2)
  temp = [ptSeries(:,iVox) params.seedTC{iScan}];
  corTemp = corr(temp); clear temp
  corVal(iVox) = corTemp(1,2); clear corTemp
end


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

