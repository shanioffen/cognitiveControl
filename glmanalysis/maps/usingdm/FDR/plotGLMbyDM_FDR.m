function plotGLMbyDM_FDR(view,overlayNum,scan,realX,realY,realS,roi) %,xBase,yBase,sBase)
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
d.pvalTH = analysis.overlays(3).data{scan};
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

% get pvalTH cutoff value
viewSet(view,'currentOverlay',3); % pvalTH is the third overlay
cutoffpvalTH = viewGet(view,'overlayMin');
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
  if(length(d.stimNames) == 4)
    stimNames{1} = 's1';
    stimNames{2} = 'correctDelay';
    stimNames{3} = 'wrongDelay';
    stimNames{3} = 's2';
  elseif length(d.stimNames)==6
    stimNames{1} = 's1correct';
    stimNames{2} = 'correctDelay';
    stimNames{3} = 's2correct';
    stimNames{4} = 's1wrong';
    stimNames{5} = 'wrongDelay';
    stimNames{6} = 's2wrong';
  else
    for istim = 1:length(d.stimNames)
      stimNames{istim} = num2str(d.stimNames(istim));
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
	
	if (d.pvalTH(x,y,s) >= cutoffpvalTH) & (overlayData(x,y,s) >= overlayCutoff) % only take voxels that have sig r-2 and meet current overlay cutoff
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
      if d.pvalTH(x,y,s) >= cutoffpvalTH
	roin(roinum) = roin(roinum)+1;
	[ehdr(roin(roinum),:,:) time] = gethdr(fixd,x,y,s);
	multiBetas(roin(roinum),:) = squeeze(fixd.ehdr(x,y,s,:));
      end
    end
    % plot the average of the ehdrs that beat the pvalTH cutoff
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
% to do this, need the experiment timing, and to do that, you need to know
% the subject and experiment. 
prompt={'Subject name','Expt name'};
name='enter subj and expt';
subjExpt=inputdlg(prompt,name);

global subj;
global expt;
subj = subjExpt{1};
expt = subjExpt{2};
exptTiming = getExptTiming(subj,expt);

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
tmp1 = loadTSeries(dataView,scan,realS,[],realX,realY);
tmp = squeeze(tmp1);
tSeries = percentTSeries(tmp);
clear tmp tmp1
disp('done')

% get trial triggered average
trialTrigSeries = getTTA(tSeries,exptTiming);
trialTrigSeries(:,3:4) = getTTA(tSeries,exptTiming,0);
% calculate model
modelSeries = DM*betas; 

% get TTA of model
trialTrigModel(:,1:2) = getTTA(modelSeries,exptTiming);

% get R-Squared
[rsq res] = calcVarAccnt(trialTrigSeries(:,1),trialTrigModel(:,1));
disp('should be keyboard now')
keyboard
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
% plot correct vs incorrect
if size(trialTrigSeries,2)==4
  figure(fignum), subplot(3,2,5);
  legendHandle(1) = errorbar(-4:2:30,trialTrigSeries(:,1),trialTrigSeries(:,2),'r.-.');
  legendStr{1} = 'Correct';
  hold on
  legendHandle(2) = errorbar(-4:2:30,trialTrigSeries(:,3),trialTrigSeries(:,4),'b.-.');
  legendStr{2} = 'Incorrect';
  xlabel('time');
  xlim([-4 30]);
  ylabel('MRI signal');
  legend(legendHandle,legendStr);
  title(sprintf('cube TTA (%i,%i,%i), Rsqu = %0.2f',xCoord,yCoord,sCoord,rsq));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the trial-triggered averages for data and model: for the ROI or cube around the voxel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check if user really wants to spend the time
ifPlot = questdlg('Do you want to plot the TTA for the ROI or cube?');
if strcmp(ifPlot,'No')||strcmp(ifPlot,'Cancel')
  return
end
clear ifPlot

% keyboard % may not want to take the time to do this

tic
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
	
	if (d.pvalTH(x,y,s) >= cutoffpvalTH) & (overlayData(x,y,s) >= overlayCutoff)
	  count = count+1;
	  % load that voxel's time course
          disp(sprintf('loading cube tseries %i',count))
          tmp1 = loadTSeries(dataView,scan,s,[],x,y);
          tmp = squeeze(tmp1);
          cubeTseries(:,count) = percentTSeries(tmp);
          clear tmp tmp1
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
    trialTrigCube(:,1:2) = getTTA(cubeTseries,exptTiming);
    trialTrigCube(:,3:4) = getTTA(cubeTseries,exptTiming,0);

    % calculate the model:
    modelCube = DM*Mbetas'; % the model for the cube is the DM * the avge betas - 
			    
    % get TTA of the model:				
    trialTrigModCube = getTTA(modelCube,exptTiming);
    
    % get R-Squared
    [rsq res] = calcVarAccnt(trialTrigCube(:,1),trialTrigModCube(:,1));
    
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
      trialTrigROI = getTTA(roiTseries,exptTiming);
      trialTrigROI(:,3:4) = getTTA(roiTseries,exptTiming,0);
      
      % calculate ROI model
      modelROI = DM*Mbetas'; % the model for the ROI is the DM * the avge betas - 
      
      % get TTA of ROI model				 
      trialTrigModROI = getTTA(modelROI,exptTiming);

      % get R-Squared
      [rsq res] = calcVarAccnt(trialTrigROI(:,1),trialTrigModROI(:,1));
      
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
      % plot correct vs incorrect
      if size(trialTrigROI,2)==4
        figure(fignum), subplot(3,2,6);
        legendHandle(1) = errorbar(-4:2:30,trialTrigROI(:,1),trialTrigROI(:,2),'r.-.');
        legendStr{1} = 'Correct';
        hold on
        legendHandle(2) = errorbar(-4:2:30,trialTrigROI(:,3),trialTrigROI(:,4),'b.-.');
        legendStr{2} = 'Incorrect';
        xlabel('time');
        xlim([-4 30]);
        ylabel('MRI signal');
        legend(legendHandle,legendStr);
        title(sprintf('cube TTA (%i,%i,%i), Rsqu = %0.2f',xCoord,yCoord,sCoord,rsq));
      end

    end
  end
end

t = toc;
disp(sprintf('took %i seconds',t))
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
% can have 3 predictors (s1,d,s2) 
% or 6 predictors (s1,d,s2 correct, s1, d, s2 wrong)
% note this only works now for when model trials separately for correct and incorrect, or all together
% ie, can't do just modeling delays separately, got too lazy to keep track of everything
if numBetas == 4, disp('Error plotting: plotGLMbyDM_so doesnt work for modeling only delays anymore'), return, end

% the timing is in the analysis structure:
trialStarts = find(d.DM(:,1));
numTrials = length(trialStarts);
epochTseries = NaN*ones(numTrials,18); % 18 time points = the longest delay period (8) + ITI  (8) + 2 volumes before trial start
cutoff = 4; % only include 3 longest delays in graph; delays will only register as 1-6 TRs because of how DM is made

% need to be able to get stim durations, which requires a little extra coding right now:
delays = d.DM(:,2);
for itrial = 1:length(trialStarts)-1
  sumstart = trialStarts(itrial);
  sumend = min(length(delays),trialStarts(itrial+1)-1); % just end before next start
  delayDuration(itrial) = sum(delays(sumstart:sumend));
end
delayDuration(end+1) = sum(delays(trialStarts(end):end));
% if model S2 as starting during delay, need to add 1 to these durations
trialEnds = find(d.DM(:,3));
if length(trialEnds) > 2* length(trialStarts)
  delayDuration = delayDuration + 1;
end
delayDuration = delayDuration';

numValidTrials = length(find(delayDuration>=cutoff));

% gather the trials with durations above the cutoff
for iTrial = 1:numTrials
  if delayDuration(iTrial)>=cutoff
    startPnt = max(1,trialStarts(iTrial)-2);  % take 2 TRs before trial starts
    endPnt = min(startPnt + 17,length(tSeries));
    epochLength = endPnt-startPnt+1;
    epochTseries(iTrial,1:epochLength) = tSeries(startPnt:endPnt);
  end
end

trialTrigSeries(:,1) = nanmean(epochTseries);
trialTrigSeries(:,2) = nanstd(epochTseries)/sqrt(numValidTrials);
keyboard
% plot also the incorrect trials if modeled separately -
if numBetas == 6
  clear delayDuration numValidTrials 
  clear startPnt endPnt epochLength epochTseries
  clear trialStarts numTrials;
  clear sumstart sumend delays
  
  trialStarts = find(d.DM(:,4));
  numTrials = length(trialStarts);
  epochTseries = NaN*ones(numTrials,18); 

  % need to be able to get stim durations, which requires a little extra coding right now:
  delays = d.DM(:,5);
  for itrial = 1:length(trialStarts)-1
    sumstart = trialStarts(itrial);
    sumend = min(length(delays),trialStarts(itrial+1)-1); % just end before next start
    delayDuration(itrial) = sum(delays(sumstart:sumend));
  end
  delayDuration(end+1) = sum(delays(trialStarts(end):end));
  clear sumstart sumend delays;
  
  % if model S2 as starting during delay, need to add 1 to these durations
  trialEnds = find(d.DM(:,6));
  if length(trialEnds) > 2* length(trialStarts)
    delayDuration = delayDuration + 1;
  end
  delayDuration= delayDuration';
  
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


