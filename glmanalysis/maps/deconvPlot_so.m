function deconvPlot_so(view,overlayNum,scan,realX,realY,realS,roi)
% eventRelatedPlot.m
%
%       $Id: eventRelatedPlot.m,v 1.23 2007/09/13 21:21:33 farshadm Exp $	
%      usage: eventRelatedPlot()
%         by: justin gardner
%       date: 10/20/06
%    purpose: 
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

% check arguments
if ~any(nargin == [1:7])
  help eventRelatedPlot
  return
end

% get the analysis structure
analysis = viewGet(view,'analysis');
d = analysis.d{scan};
if isempty(d)
  disp(sprintf('(eventRelatedPlot) Event related not for scan %i',scan));
  return
end
d.r2 = analysis.overlays(1).data{scan};
% select the window to plot into
selectGraphWin;

global MLR;
fignum = MLR.graphFigure+1;

% turn off menu/title etc.
%set(fignum,'NumberTitle','off');
%set(fignum,'Name','eventRelatedPlot');

% set roi coords
for roinum = 1:length(roi)
  % get scan coordinates
  roi{end}.scanCoords = getROICoordinates(view,roi{roinum},scan);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ehdr time ehdrste] = gethdr(d,realX,realY,realS);
figure(fignum)
subplot(2,2,1);
% display ehdr with out lines if we have a fit
% since we also need to plot fit
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
title(sprintf('Voxel (%i,%i,%i): r2=%0.3f',realX,realY,realS,analysis.overlays(1).data{scan}(realX,realY,realS)));
xaxis(0,max(time));
% add peaks if they exist to the legend
if isfield(d,'stimNames')
  stimNames = d.stimNames;
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
  ehdr = [];
  count = 0;
  maxVals = MLR.groups(groupNum).scanParams.dataSize;
  for xRange = -1:1
    x = realX + xRange;
    if x<1 | x>maxVals(1), continue, end 

    for yRange = -1:1
      y = realY + yRange;
      if y<1 | y>maxVals(2), continue, end 

      for sRange = -1:1
	s = realS + sRange;
	if s<1 | s>maxVals(3), continue, end 
	
	if (d.r2(x,y,s) >= cutoffr2) & (overlayData(x,y,s) >= overlayCutoff) 
        % only take voxels that have sig r-2 and meet current overlay cutoff
        % this also restricts the cube to fxnal voxels
	  count = count+1;
	  [ehdr(count,:,:) time] = gethdr(d,x,y,s);
	end
      end
    end
  end
  
  % plot the average of the ehdrs that beat the r2 cutoff
  if count
    figure(fignum), subplot(2,2,2);
    plotEhdr(time,squeeze(mean(ehdr)),squeeze(std(ehdr))/sqrt(count));
    title(sprintf('Cube (n=%i/%i)',count,27));
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
	[ehdr(roin(roinum),:,:) time] = gethdr(d,x,y,s);
      end
    end
    % plot the average of the ehdrs that beat the r2 cutoff
    if roin(roinum)
      figure(fignum), subplot(length(roi)+1,2,2);
      plotEhdr(time,squeeze(mean(ehdr)),squeeze(std(ehdr))/sqrt(size(roi{roinum}.scanCoords,2)));
      DPI = multiBetas(:,2)./multiBetas(:,1);
      title(sprintf('%s (n=%i/%i)',roi{roinum}.name,roin));
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


