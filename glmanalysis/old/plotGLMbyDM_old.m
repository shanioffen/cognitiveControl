function plotGLMbyDM_so(view,overlayNum,scan,realX,realY,realS,roi) %,xBase,yBase,sBase)
% eventRelatedPlot.m
%
%      usage: eventRelatedPlot_so()
%         by: shani, adapted from code by justin gardner
%       date: August 2007
%    purpose: plots the betas and the trial-triggered average for a voxel and surrounding voxels (or surrounding ROI)
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
d = analysis.d{scan};
d.r2 = analysis.overlays(1).data{scan};
DM = d.scm;

% select the window to plot into
selectGraphWin;

global MLR;
fignum = MLR.graphFigure;

% turn off menu/title etc.
set(fignum,'NumberTitle','off');
set(fignum,'Name','eventRelatedPlot');

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

% get R2 cutoff value (since current overlay may not be R2
viewSet(view,'currentOverlay',1); % r2 is the first overlay
cutoffr2 = viewGet(view,'overlayMin');
viewSet(view,'currentOverlay',currOverNum); % set it back to whatever the current overlay really is

if isempty(d)
  disp('No analysis');
  return
end

% set the group number
groupNum = viewGet(view,'currentGroup');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the hemodynamic response for voxel
% this is really the canonical hrf * the betas for each component
% ** This code also gets the betas which are used 
% later in the model graphing code **
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if running all 5 subjects need to choose one hrf to display
fixd = d; fixd.hrf = d.hrf(:,1);

figure(fignum), subplot(max(length(roi)+2,3),2,1);
[ehdr time ehdrste] = gethdr(fixd,realX,realY,realS);
betas = squeeze(d.ehdr(realX,realY,realS,:));

% I don't know what this peak and fit stuff is
if isfield(d,'peak') & isfield(d.peak,'fit') & ~any(isnan(d.peak.amp(realX,realY,realS,:)))
  plotEhdr(time,ehdr,ehdrste,'');
  for r = 1:d.nhdr
    d.peak.fit{realX,realY,realS,r}.smoothX = 1:.1:d.hdrlen;
    fitTime = d.tr*(d.peak.fit{realX,realY,realS,r}.smoothX-0.5);
    plot(fitTime+d.tr/2,d.peak.fit{realX,realY,realS,r}.smoothFit,getcolor(r,'-'));
  end
else
  plotEhdr(time,ehdr,ehdrste);
end
title(sprintf('Voxel (%i,%i,%i): r2=%0.2f, betas = %0.1f, %0.1f, DPI = %0.2f',...
	      xCoord,yCoord,sCoord,analysis.overlays(1).data{scan}(realX,realY,realS),betas(1),betas(2), betas(2)/betas(1)));
xaxis(0,max(time));
% add peaks if they exist to the legend
if isfield(d,'stimNames')
  stimNames = num2str(d.stimNames);
  if isfield(d,'peak')
    for i = 1:d.nhdr
      stimNames{i} = sprintf('%s: %s=%0.2f',stimNames{i},d.peak.params.method,d.peak.amp(realX,realY,realS,i));
    end
  end
  legend(stimNames);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if there is an roi at this voxel
% then plot mean (beta*HRF)
%
% If not, then plot the mean response
% of a cube of voxels centered on this voxel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(length(roi)==0)
  figure(fignum), subplot(3,2,2);
  ehdr = [];
  count = 0;
  maxVals = MLR.groups(groupNum).scanParams.dataSize;
  for xRange = -1:1
    x = realX + xRange;
    if x<1 | x>maxVals(1), continue, end % make sure we haven't gone out of the range of the data

    for yRange = -1:1
      y = realY + yRange;
      if y<1 | y>maxVals(2), continue, end % make sure we haven't gone out of the range of the data

      for sRange = -1:1
	s = realS + sRange;
	if s<1 | s>maxVals(3), continue, end % make sure we haven't gone out of the range of the data
	
	if (d.r2(x,y,s) >= cutoffr2) & (overlayData(x,y,s) >= overlayCutoff) % only take voxels that have sig r-2 and meet current overlay cutoff
									     % this also restricts the cube to fxnal voxels
	  count = count+1;
	  [ehdr(count,:,:) time] = gethdr(fixd,x,y,s);
	  multiBetas(count,:) = squeeze(d.ehdr(x,y,s,:));
	end
      
      end
    end
  end
  
  % plot the average of the ehdrs that beat the r2 cutoff
  if count
    plotEhdr(time,squeeze(mean(ehdr)),squeeze(std(ehdr))/sqrt(count));
    Mbetas = mean(multiBetas,1);
    DPI = multiBetas(:,2)./multiBetas(:,1);
    
    title(sprintf('Cube (n=%i/%i), mean betas= %0.1f %0.1f, mean DPI = %0.2f',...
		  count,27, Mbetas(1), Mbetas(2), Mbetas(2)/Mbetas(1)), 'Interpreter','none');
  end

else % if there are ROIs
  
  for roinum = 1:length(roi)
    ehdr = [];
    roin(roinum) = 0;
    for voxnum = 1:size(roi{roinum}.scanCoords,2)
      x = roi{roinum}.scanCoords(1,voxnum);
      y = roi{roinum}.scanCoords(2,voxnum);
      s = roi{roinum}.scanCoords(3,voxnum);
      if d.r2(x,y,s) >= cutoffr2
	roin(roinum) = roin(roinum)+1;
	[ehdr(roin(roinum),:,:) time] = gethdr(fixd,x,y,s);
	multiBetas(roin(roinum),:) = squeeze(fixd.ehdr(x,y,s,:));
      end
    end
    % plot the average of the ehdrs that beat the r2 cutoff
    if roin(roinum)
      figure(fignum), subplot(length(roi)+2,2,2);
      plotEhdr(time,squeeze(mean(ehdr)),squeeze(std(ehdr))/sqrt(size(roi{roinum}.scanCoords,2)));
      Mbetas = mean(multiBetas,1);
      DPI = multiBetas(:,2)./multiBetas(:,1);
      title(sprintf('%s (n=%i/%i), mean betas= %0.2f %0.2f, mean DPI = %0.2f', ...
                    roi{roinum}.name,roin,size(roi{roinum}.scanCoords,2), Mbetas(1), Mbetas(2), Mbetas(2)/Mbetas(1)), 'Interpreter','none');
    end
  end
end
drawnow;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%
% plot the trial-triggered averages for data and model: For the chosen voxel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%

% check if user really wants to spend the time
ifPlot = questdlg('Do you want to plot the TTA?');
if strcmp(ifPlot,'No')||strcmp(ifPlot,'Cancel')
  return
end
clear ifPlot

% make sure you're loading the unblurred data
dataView = newView('Volume');
concatNum = viewGet(dataView,'groupnum','Concatenation');
if ~isempty(concatNum)
  dataView = viewSet(dataView,'currentgroup',concatNum);
else
  dataView = view;
end
%actually, I'm not sure about this, so skip it for now....
dataView = view;
  
% load tSeries at chosen voxel
disp('loading voxel t series')
tmp = squeeze(loadTSeries(dataView,scan,realS,[],realX,realY));
tSeries = percentTSeries(tmp);
clear tmp
disp('done')

% get trial triggered average
trialTrigSeries = getTrialTrigAvge(tSeries,dataView,d);

% calculate model
modelSeries = DM*betas; 

% get TTA of model
trialTrigModel = getTrialTrigAvge(modelSeries,dataView,d); 

% get R-Squared
[rsq res] = calcVarAccnt(trialTrigSeries,trialTrigModel);

% plot
figure(fignum), subplot(max(length(roi)+2,3),2,3)
legendHandle(1) = errorbar(-4:2:30,trialTrigSeries(:,1),trialTrigSeries(:,2),'r.-.');
legendStr{1} = 'Data';
hold on
legendHandle(2) = errorbar(-4:2:30,trialTrigModel(:,1),trialTrigModel(:,2),'k.-');
legendStr{2} = 'Model';
xlabel('time');
xlim([-4 30]);
ylabel('MRI signal');
title(sprintf('TTA (%i,%i,%i), rSqu = %0.2f',xCoord,yCoord,sCoord,rsq));
legend(legendHandle,legendStr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the trial-triggered averages for data and model: for the ROI or cube around the voxel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% keyboard % may not want to take the time to do this

if(length(roi)==0)

  count = 0;
  maxVals = MLR.groups(groupNum).scanParams.dataSize;
  for xRange = -1:1
    x = realX + xRange;
    if x<1 | x>maxVals(1), continue, end % make sure we haven't gone out of the range of the data

    for yRange = -1:1
      y = realY + yRange;
      if y<1 | y>maxVals(2), continue, end % make sure we haven't gone out of the range of the data

      for sRange = -1:1
	s = realS + sRange;
	if s<1 | s>maxVals(3), continue, end % make sure we haven't gone out of the range of the data
	
	if (d.r2(x,y,s) >= cutoffr2) & (overlayData(x,y,s) >= overlayCutoff)
	  count = count+1;
	  % load that voxel's time course
          disp(sprintf('loading cube tseries %i',count))
          tmp = squeeze(loadTSeries(dataView,scan,s,[],x,y));
          cubeTseries(:,count) = percentTSeries(tmp);
          clear tmp
          disp('done')

	end
      end
    end
  end  
  
  % If anything in the cube exceeds the cutoff:
  if(count)
    % get avge time course for the cube:
    cubeTseries = mean(cubeTseries,2);
    
    % get trial triggered average for this avge time course:
    trialTrigCube = getTrialTrigAvge(cubeTseries,dataView,d);

    % calculate the model:
    modelCube = DM*Mbetas'; % the model for the cube is the DM * the avge betas - 
			    
    % get TTA of the model:				
    trialTrigModCube = getTrialTrigAvge(modelCube,dataView,d);  
    
    % get R-Squared
    [rsq res] = calcVarAccnt(trialTrigCube,trialTrigModCube);
    
    % plot:
    figure(fignum), subplot(3,2,4);
    legendHandle(1) = errorbar(-4:2:30,trialTrigCube(:,1),trialTrigCube(:,2),'r.-.');
    legendStr{1} = 'Data';
    hold on
    legendHandle(2) = errorbar(-4:2:30,trialTrigModCube(:,1),trialTrigModCube(:,2),'k.-');
    legendStr{2} = 'Model';
    xlabel('time');
    xlim([-4 30]);
    ylabel('MRI signal');
    legend(legendHandle,legendStr);
    title(sprintf('cube TTA (%i,%i,%i), Rsqu = %0.2f',xCoord,yCoord,sCoord,rsq))
  
    % plot correct vs incorrect
    if size(trialTrigCube,2)==4
      figure(fignum), subplot(3,2,6);
      legendHandle(1) = errorbar(-4:2:30,trialTrigCube(:,1),trialTrigCube(:,2),'r.-.');
      legendStr{1} = 'Correct';
      hold on
      legendHandle(2) = errorbar(-4:2:30,trialTrigCube(:,3),trialTrigCube(:,4),'b.-.');
      legendStr{2} = 'Incorrect';
      xlabel('time');
      xlim([-4 30]);
      ylabel('MRI signal');
      legend(legendHandle,legendStr);
      title(sprintf('cube TTA (%i,%i,%i), Rsqu = %0.2f',xCoord,yCoord,sCoord,rsq));
    end
  end

else % if there are ROIs
  
  for roinum = 1:length(roi)
    ehdr = [];
    % The ROI should already be restricted to relevant voxels, don't do the restriction here, but make sure there are voxels to include
    if roin(roinum)
      
      roiTseries = tseriesROI(dataView, groupNum, roinum, scan);
      roiTseries = mean(roiTseries{1},2);
      
      % get TTA of ROI tSeries
      trialTrigROI = getTrialTrigAvge(roiTseries,dataView,d);
      
      % calculate ROI model
      modelROI = DM*Mbetas'; % the model for the ROI is the DM * the avge betas - 
                             % eventually will be better to actually do the GLM on the ROI instead of taking avge betas ********* TODO
      
      % get TTA of ROI model				 
      trialTrigModROI = getTrialTrigAvge(modelROI,dataView,d);
      
      % get R-Squared
      [rsq res] = calcVarAccnt(trialTrigROI,trialTrigModROI);
      
      % plot
      figure(fignum), subplot(length(roi)+2,2,(roinum+3));
      legendHandle(1) = errorbar(-4:2:30,trialTrigROI(:,1),trialTrigROI(:,2),'r.-.');
      legendStr{1} = 'Data';
      hold on
      legendHandle(2) = errorbar(-4:2:30,trialTrigModROI(:,1),trialTrigModROI(:,2),'k.-');
      legendStr{2} = 'Model';
      xlabel('time');
      xlim([-4 30]);
      ylabel('MRI signal');
      legend(legendHandle,legendStr);
      title(sprintf('TTA %s, Rsqu = %0.2f',roi{roinum}.name,rsq));

    end
  end
end
drawnow;


%%%%%%%%%%%%%%%%%%%%%%%%%
% function to plot ehdr
%%%%%%%%%%%%%%%%%%%%%%%%%
function plotEhdr(time,ehdr,ehdrste,lineSymbol)

% whether to plot the line inbetween points or not
if ~exist('lineSymbol','var'),lineSymbol = '-';,end

% and display ehdr
for i = 1:size(ehdr,1)
  if nargin == 2
    h=plot(time,ehdr(i,:),getcolor(i,getsymbol(i,lineSymbol)),'MarkerSize',8);
  else
    h=errorbar(time,ehdr(i,:),ehdrste(i,:),ehdrste(i,:),getcolor(i,getsymbol(i,lineSymbol)),'MarkerSize',8);
  end
  set(h,'MarkerFaceColor',getcolor(i));
  hold on
end
xlabel('Time (sec)');
ylabel('% Signal change');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to get Trial Triggered Average
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function trialTrigSeries = getTrialTrigAvge(tSeries,view,d);

% need to know the model
numBetas = size(d.DM,2);
% can have 3 predictors (s1,d,s2) or 4 predictors (s1, delay_correct, delay_wrong, s2) 
% or 6 predictors (s1,d,s2 correct, s1, d, s2 wrong)

% the timing is in the analysis structure:
trialStarts = find(d.DM(:,1));
numTrials = length(trialStarts);
epochTseries = NaN*ones(numTrials,18); % 18 time points = the longest delay period (8) + ITI  (8) + 2 volumes before trial start
cutoff = 6; % delay must be more than 12 sec to go in graph

% need to be able to get stim durations, which requires a little extra coding right now:
% if the scan ended in the middle of the delay period, the first and last stim components
% won't be aligned, so need to fix that:
if numBetas == 3 || numBetas == 6 % s1 d s2 or s1correct dcorrect s2correct
  s1place = 1;
  dplace = 2;
  s2place = 3;
elseif numBetas == 4 % s1 dcorrect dwrong s2
  s1place = 1;
  dplace = 2;
  s2place = 4;
end

if sum(d.DM(:,s2place)) > sum(d.DM(:,s1place))
  tempEnds = find(d.DM(:,s2place));
  trialEnds = tempEnds(1:2:length(tempEnds));
else
  trialEnds = find(d.DM(:,s2place));
end


numSkips = length(trialStarts) - length(trialEnds);
comp1 = trialStarts;
comp3 = NaN*ones(size(comp1));
comp3(1:length(trialEnds)) = trialEnds;

if numSkips > 0 
  for iSkip = 1:numSkips
    skipInd = find(comp3 - comp1 > 7);
    comp3(skipInd(1)+1:end) = comp3(skipInd(1):end-1);
    comp3(skipInd(1)) = NaN;
  end
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

% plot also the incorrect trials if modeled separately -
clear numSkips comp1 comp3 skipInd delayDuration numValidTrials 
clear startPnt endPnt epochLength epochTseries
clear trialStarts numTrials s1times s2times

trialStarts = find(d.DM(:,4));
numTrials = length(trialStarts);
epochTseries = NaN*ones(numTrials,18); 

if numBetas == 6
  if sum(d.DM(:,6)) > sum(d.DM(:,4))
    tempEnds = find(d.DM(:,6));
    trialEnds = tempEnds(1:2:length(tempEnds));
  else
    trialEnds = find(d.DM(:,6));
  end
  
  numSkips = length(trialStarts) - length(trialEnds);
  comp1 = trialStarts;
  comp3 = NaN*ones(size(comp1));
  comp3(1:length(trialEnds)) = trialEnds;
  if numSkips > 0 
    for iSkip = 1:numSkips
      skipInd = find(comp3 - comp1 > 7);
      comp3(skipInd(1)+1:end) = comp3(skipInd(1):end-1);
      comp3(skipInd(1)) = NaN;
    end
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
  
  trialTrigSeries(:,3) = nanmean(epochTseries);
  trialTrigSeries(:,4) = nanstd(epochTseries)/sqrt(numValidTrials);

end

