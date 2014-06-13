%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TTA = getTTA(TC,timing,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the TTA for a time course, given the timing of the expt
% can get TTA only for correct trials (right = 1), only incorrect (right = -1), 
% or all (right = 0) Default is 1
  
if nargin < 3
  right = 1;
  numBins = 1;
elseif nargin == 3
  right = varargin{1};
  numBins = 1;
elseif nargin == 4
  right = varargin{1};
  numBins = varargin{2};
end

% find which trials meet criteria:
correct = timing.correct;
if right == 1
  takeIndx = find(correct == 1);
elseif right == -1
  takeIndx = find(correct == 0);
elseif right == 0
  takeIndx = 1:length(correct);
end
  

% run through the bins
for iBin = 1:numBins
  if numBins ==1
    upperLim = 8;
    lowerLim = 5;
  else
    binSize = 8/numBins;
    binList = binSize:binSize:8;
    upperLim = binList(iBin);
    lowerLim = upperLim - binSize;
  end % check if only 1 bin
  
  % check the delays for those trials
  starts = timing.startTR(takeIndx);
  delays = timing.delayTR(takeIndx);

  keepIndx = find(and(delays <= upperLim, delays > lowerLim));
  starts = starts(keepIndx);
  
  numTrials = length(keepIndx);
  epochTseries = NaN*ones(numTrials,18); 
  % 18 time points = the longest delay period (8) + ITI  (8) + 2 volumes before trial start
  
  takePnts = 2 + 8 + upperLim - 1;
  for iTrial = 1:numTrials
    if starts(iTrial) == 1 % for first trial, special case
      startPnt = 1;
      endPnt = startPnt + takePnts - 2; % since not get 2 pnts before trial start
      epochTseries(iTrial,3:takePnts+1) = TC(startPnt:endPnt);
    else
      startPnt = starts(iTrial)-2; 
      endPnt = min(startPnt + takePnts, length(TC)); % don't go off other edge either 
      epochLength = endPnt - startPnt + 1;
      epochTseries(iTrial,1:epochLength) = TC(startPnt:endPnt);
    end
  end
  
  % there's a problem with how I'm doing this: when getting the final trial before a scan
  % ended, I'm going to be taking data starting in the next scan too when I take all 18 pnts
  % might need to go through things scan by scan, but that's a pain. Or figure out a way
  % to easily calculate where the transitions happen, by subtracting starts from itself
  % shifted one, and then deal with that.
  TTA(:,2*(iBin-1)+1) = nanmean(epochTseries);
  TTA(:,2*iBin) = nanstd(epochTseries)./sqrt(numTrials);
  
  clear starts delays keepIndx numTrials epochTseries takePnts

end % for iBin
