function corPlot_so(view,overlayNum,scan,realX,realY,realS,roi) %,xBase,yBase,sBase)
% corPlot_so.m
%
%      usage: eventRelatedPlot_so()
%         by: shani, adapted from code by justin gardner
%       date: October 2007
%    purpose: plots the the trial-triggered average for a voxel and surrounding voxels (or surrounding ROI)
%             plots the TTA for the seed as well, and gives the corr to of the voxel to the seed
%


baseCoords = viewGet(view,'mouseDownBaseCoords');
xBase = baseCoords(1); yBase = baseCoords(2); sBase = baseCoords(3);
base2tal = viewGet(view,'base2tal'); % keyboard
if(~isempty(base2tal)) % if there is a tal transform, use Tal coordinates in titles
  talCoords = round(base2tal * [xBase yBase sBase 1]');
  xTal = talCoords(1); yTal = talCoords(2); zTal = talCoords(3);
  xCoord = xTal; yCoord = yTal; sCoord = zTal;
else
  xCoord = xBase; yCoord = yBase; sCoord = sBase;
end

% get the analysis structure
analysis = viewGet(view,'analysis');
seedTC = analysis.params.seedTC;
seedType = analysis.params.seedType;
d = analysis.d{scan};
corVal = analysis.overlays(1).data{scan}(realX,realY,realS);
allCorVals = analysis.overlays(1).data{scan};

% select the window to plot into
selectGraphWin;

global MLR;
fignum = MLR.graphFigure;

% turn off menu/title etc.
set(fignum,'NumberTitle','off');
set(fignum,'Name','CorValPlot');

% set roi coords
for roinum = 1:length(roi)
  % get scan coordinates
  roi{roinum}.scanCoords = getROICoordinates(view,roi{roinum},scan);
end

% get current cutoff value
currOverName = viewGet(view,'overlayName');
currOverNum = viewGet(view,'currentOverlay'); % save for later
overlayCutoff = viewGet(view,'overlayMin');
overlayData = analysis.overlays(currOverNum).data{scan};

if isempty(d)
  disp('No analysis');
  return
end

% set the group number
groupNum = viewGet(view,'currentGroup');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%
% plot the trial-triggered averages for data and model: For the chosen voxel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%

figure(fignum), subplot(max(length(roi),1),2,1)

% load tSeries at chosen voxel
tSeries = percentTSeries(squeeze(loadTSeries(view,scan,realS,[],realX,realY)));

% get trial triggered average
trialTrigSeries = getTrialTrigAvge(tSeries,view,d);

% get TTA of seed
trialTrigSeed = getTrialTrigAvge(seedTC{1},view,d);  % need to fix so access the right scan in the seed

% get R-Squared for seed fit to time course
[rsq res] = calcVarAccnt(trialTrigSeries,trialTrigSeed);

% plot
legendHandle(1) = errorbar(-2:2:32,trialTrigSeries(:,1),trialTrigSeries(:,2),'r.-.');
legendStr{1} = 'Chosen Vox';
hold on
legendHandle(2) = errorbar(-2:2:32,trialTrigSeed(:,1),trialTrigSeed(:,2),'k.-');
legendStr{2} = ['seed ' seedType];
xlabel('Time (s)');
ylabel('MRI signal');
title(sprintf('TTA: Vox (%i,%i,%i), corVal = %0.2f, rSqu = %0.2f',xCoord,yCoord,sCoord,corVal,rsq));
legend(legendHandle,legendStr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the trial-triggered averages for data and model: for the ROI or cube around the voxel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(length(roi)==0) % if no ROIs, plot for a cube of voxels around the chosen voxel

  figure(fignum), subplot(1,2,2);
  count = 0;
  maxVals = MLR.groups(groupNum).scanParams.dataSize;
  for xRange = -1:1 % -2:2
    x = realX + xRange;
    if x<1 | x>maxVals(1), continue, end % make sure we haven't gone out of the range of the data

    for yRange = -1:1 % -2:2
      y = realY + yRange;
      if y<1 | y>maxVals(2), continue, end % make sure we haven't gone out of the range of the data

      for sRange = -1:1 % -2:2
	s = realS + sRange;
	if s<1 | s>maxVals(3), continue, end % make sure we haven't gone out of the range of the data
	
	if (overlayData(x,y,s) >= overlayCutoff)
	  count = count+1;
	  % load that voxel's time course
	  cubeTseries(:,count) = percentTSeries(squeeze(loadTSeries(view,scan,s,[],x,y)));
	  
	  % calculate avge corVal
	  multiCorVal(count,1) = allCorVals(x,y,s);
	
	end
      end
    end
  end  

  if(count) % if anything in the cube exceeds the cutoff
  
    % get avge corVal for the cube:
    avgCorVal = mean(multiCorVal);
  
    % get avge time course for the cube:
    cubeTseries = mean(cubeTseries,2);
  
    % get trial triggered average for this avge time course:
    trialTrigCube = getTrialTrigAvge(cubeTseries,view,d);
    
    % get R-Squared for Cube fit to seed
    [rsq res] = calcVarAccnt(trialTrigCube,trialTrigSeed);
    
    % plot:
    legendHandle(1) = errorbar(-2:2:32,trialTrigCube(:,1),trialTrigCube(:,2),'r.-.');
    legendStr{1} = 'Chosen Vox';
    hold on
    legendHandle(2) = errorbar(-2:2:32,trialTrigSeed(:,1),trialTrigSeed(:,2),'k.-');
    legendStr{2} = ['seed ' seedType];
    xlabel('Time (s)');
    ylabel('MRI signal');
    legend(legendHandle,legendStr);
    
    title(sprintf('TTA: center %i:%i:%i (n=%i/%i), corVal = %0.2f, Rsqu = %0.2f',...
                  xCoord, yCoord, sCoord, count, 27, avgCorVal, rsq))
  end
  
else % if there are ROIs
  
  for roinum = 1:length(roi)
    figure(fignum), subplot(length(roi),2,(roinum+1));
    ehdr = [];
    % The ROI should already be restricted to relevant voxels, don't do the restriction here
    roiTseries = tseriesROI(view, groupNum, roinum, scan);
    roiTseries = mean(roiTseries{1},2);
    
    % get TTA of ROI tSeries
    trialTrigROI = getTrialTrigAvge(roiTseries,view,d);

    % get R-Squared of ROI fit to seed
    [rsq res] = calcVarAccnt(trialTrigROI,trialTrigSeed);

    % plot
    legendHandle(1) = errorbar(-2:2:32,trialTrigROI(:,1),trialTrigROI(:,2),'r.-.');
    legendStr{1} = 'Data';
    hold on
    legendHandle(2) = errorbar(-2:2:32,trialTrigSeed(:,1),trialTrigSeed(:,2),'k.-');
    legendStr{2} = seedType;
    xlabel('Time (s)');
    ylabel('MRI signal');
    legend(legendHandle,legendStr);

    title(sprintf('Trial-triggered average for %s, Rsqu = %0.2f',roi{roinum}.name,rsq));
    
  end
end
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to get Trial Triggered Average
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function trialTrigSeries = getTrialTrigAvge(tSeries,view,d,cutoff);

if nargin == 3
  cutoff = 6;
end


% the timing is in the analysis structure:
trialStarts = d.stimvol{1};
numTrials = length(d.stimvol{1});
epochTseries = NaN*ones(numTrials,18); % 18 time points = the longest delay period + ITI + 2 volumes before trial start

% need to be able to get stim durations, which requires a little extra coding right now:
% if the scan ended in the middle of the delay period, the first and third components
% won't be aligned, so need to fix that:

numSkips = length(d.stimvol{1})-length(d.stimvol{3}); 
comp1 = d.stimvol{1};
comp3 = NaN*ones(size(comp1));
comp3(1:length(d.stimvol{3})) = d.stimvol{3};
for iSkip = 1:numSkips
  skipInd = find(comp3 - comp1 > 7);
  comp3(skipInd(1)+1:end) = comp3(skipInd(1):end-1);
  comp3(skipInd(1)) = NaN;
end
delayDuration = comp3 - comp1;
numValidTrials = length(find(delayDuration>=cutoff));

% gather the trials with durations above the cutoff
for iTrial = 1:numTrials
  if delayDuration(iTrial)>=cutoff
    startPnt = max(1,trialStarts(iTrial)-2);  
    endPnt = min(startPnt + 17,length(tSeries));
    epochLength = endPnt-startPnt+1;
    epochTseries(iTrial,1:epochLength) = tSeries(startPnt:endPnt);
  end
end

trialTrigSeries(:,1) = nanmean(epochTseries);
trialTrigSeries(:,2) = nanstd(epochTseries)/sqrt(numValidTrials);

