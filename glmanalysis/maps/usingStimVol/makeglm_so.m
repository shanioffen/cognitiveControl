% makeglm.m
%
%      usage: makeglm(d,hrf,subjRunList)
%         by: farshad moradi
%       date: 06/14/07
%       e.g.: makeglm(d, hrf)
%    purpose: makes a stimulation convolution matrix
%             for data series. must have getstimtimes already
%             run on it. if d.hrf is set then you can
%             just pass in the data structure
%
% 2008March19 - I edited FM's code so that can deal with fixed effects,
% by having different subjects with different HRFs on different runs
% subjRunList is an optional argument that tells you how many runs each subj
% did, so the correct HRF can be applied for each
%
function d = makeglm_so(d,hrf,subjRunList)

if (nargin == 1)
  if (isfield(d,'hrf'))
      hrf = d.hrf;
  else
      help makeglm;
      return
  end
elseif (nargin == 2)
  subjRunList = [];
elseif nargin == 3
  if length(subjRunList)
    if size(hrf,2) ~= length(subjRunList)
      disp('Mismatch: hrf and subjRunList dont have same number of subjects')
      return
    end
  end
else
  help makeglm;
  return
end

% make sure hrf starts from zero
hrf = hrf-repmat(hrf(1,:), size(hrf,1), 1);
% if we have only a single run then we set
% the runTransitions for that single run
if ~isfield(d,'concatInfo') || isempty(d.concatInfo)
  runTransition = [1 d.dim(4)];
  hipassfilter = [];
else
  runTransition = d.concatInfo.runTransition;
  hipassfilter = d.concatInfo.hipassfilter;
end

% go through each run of the experiment
allscm = [];
% check that runcounts match
if length(subjRunList)
  if sum(subjRunList) ~= size(runTransition,1) 
    disp('number of runs dont add up')
    return
  end
end


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
  % make stimulus convolution matrix
  for stimnum = 1:length(d.stimvol)
    % make an array containing the stimulus times
    stimarray = zeros(1,(runTransition(runnum,2)-runTransition(runnum,1)+1)*d.supersampling);
    % only use stimvols that are within this runs volume numbers
    stimarray(d.stimvol{stimnum}(find((d.stimvol{stimnum}>=runTransition(runnum,1)*d.supersampling) & ...
        (d.stimvol{stimnum}<=runTransition(runnum,2)*d.supersampling)))-runTransition(runnum,1)*d.supersampling+1) = 1/d.supersampling;
    m = convn(stimarray', hrf(:,whichSub));
    m = m(1:length(stimarray),:);
    % remove mean 
    m = m-repmat(mean(m), size(m,1), 1);
    % downsample
    m = downsample(m, d.supersampling);
    % apply the same filter as original data
    if ~isempty(hipassfilter)
        m = real(ifft(fft(m) .* repmat(hipassfilter{runnum}', 1, size(m,2)) ));
    end
    % stack stimcmatrices horizontally
    scm = [scm, m];
  end
  % stack this run's stimcmatrix on to the last one
  allscm = [allscm;scm];
end

% set values
d.nhdr = length(d.stimvol);
d.scm = allscm;
d.hdrlen = size(hrf,2); if length(subjRunList), d.hdrlen = 1; end %otherwise causes problems
d.volumes = 1:d.dim(4);
d.simulatedhrf = downsample(hrf, d.supersampling);

function ds = downsample(s, factor)
  % downsample a signal by factor. the original and downsampled signals
  % have equal integrals.
  %
  % to downsamlple conserving the amplitude use: 
  %     ds = downsample(s, factor) / factor;
  
  a = cumsum(s, 1);
  ds = diff([zeros(1, size(s,2));a(factor:factor:end,:)]);
