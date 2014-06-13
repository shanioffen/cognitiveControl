% talairach.m
%
%        $Id: talairach.m,v 1.2 2007/12/20 16:46:49 justin Exp $
%      usage: talairach(volume)
%         by: justin gardner
%       date: 05/04/07
%    purpose: 
%             program to talairach a volume, call with filename
%             talairach('jg041001');
%
function retval = talairach_jg(event,isEvent)

% check arguments
if ~any(nargin == [1 2])
    help talairach_jg
    retval = [];
    return
end

global gTalairach;

% init arguments
if nargin == 1
  % if we are passed in a structure then this is a button callback
  disp('in the right code')
  if isstruct(event)
    event = event.event;
    retval = [];
    % otherwise it is init event
  elseif isstr(event)
    filename = sprintf('%s.img',stripext(event));
    event = 'init';
    % check for file
    if isfile(filename)
      % read header
      hdr = cbiReadNiftiHeader(filename);
      % if this is a 4D volume then we have to take a particular volume or
      % or the mean. Ask the user what to do.
      doMean = 0;subset = [];
      if hdr.dim(5) > 1
	paramsInfo = {{'volNum',0,'incdec=[-1 1]','round=1',sprintf('minmax=[0 %i]',hdr.dim(5)),'Choose volume number to display (0 for mean)'}};
	params = mrParamsDialog(paramsInfo,'Choose volume to use (0 for mean)');
	if isempty(params),return,end
	if params.volNum ~= 0
	  subset = {[],[],[],params.volNum};
	else
	  doMean = 1;
	end
      end
      % read it
      disppercent(-inf,sprintf('Loading %s',filename));
      if ~isempty(subset)
	[vol hdr] = cbiReadNifti(filename,subset);
      else
	[vol hdr] = cbiReadNifti(filename);
      end
      if doMean,vol = mean(vol,4);end
      disppercent(inf);
    else
      disp(sprintf('(talairach_jg) Could not open file %s',filename));
      return
    end
  end
end

switch (event)
  case 'init'
    initHandler(filename,vol,hdr);
  case 'end'
    endHandler;
  case 'mouseMove'
     mouseMoveHandler;
  case 'mouseUp'
    mouseUpHandler;
 case 'mouseDown'
  disp('called mouseDown')
    mouseDownHandler;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mousedown
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mouseDownHandler

global gTalairach;
disp('in mouseDownHandler')
% get the pointer location
pointerLoc = get(gca,'CurrentPoint');
% get mouse, remembering that we have swapped y/x
% (see getImageSlice) 
mouseX = round(pointerLoc(1,1));
mouseY = round(pointerLoc(1,2));

% which figure we are on
if gcf == gTalairach.fig(1)
    a = [1 2 3];
    %  mouseX = size(gTalairach.vol,2)-mouseX+1;
    mouseY = size(gTalairach.vol,3)-mouseY+1;
    x = gTalairach.pos(1);
    y = mouseX;
    z = mouseY;
elseif gcf == gTalairach.fig(2)
    a = [2 1 3];
    %  mouseX = size(gTalairach.vol,1)-mouseX+1;
    mouseY = size(gTalairach.vol,3)-mouseY+1;
    x = mouseX;
    y = gTalairach.pos(2);
    z = mouseY;
elseif gcf == gTalairach.fig(3)
    a = [3 1 2];
    %  mouseX = size(gTalairach.vol,1)-mouseX+1;
    mouseY = size(gTalairach.vol,2)-mouseY+1;
    x = mouseX;
    y = mouseY;
    z = gTalairach.pos(3);
else
    return
end

% display the point being clicked on
disp(sprintf('x:%i y:%i z:%i',x,y,z));

% display the correct volume
if ((mouseX > 0) && (mouseX < gTalairach.dim(a(2))) && ...
    (mouseY > 0) && (mouseY < gTalairach.dim(a(3))))
    figure(gTalairach.fig(a(2)));
    dispVolumeSlice(a(2),mouseX);
    gTalairach.pos(a(2)) = mouseX;

    figure(gTalairach.fig(a(3)));
    dispVolumeSlice(a(3),mouseY);
    gTalairach.pos(a(3)) = mouseY;
    set(gTalairach.fig(a(1)),'pointer','fullcrosshair');
else
    set(gTalairach.fig(a(1)),'pointer','arrow');
end

if(1)% mglGetKeys(57)
  setTalairachPoint(x,y,z);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   setTalairachPoint   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setTalairachPoint(x,y,z)

global gTalairach;

paramsInfo = {};
paramsInfo{end+1} = {'whichPoint',putOnTopOfList('None',gTalairach.refPoints),'Use this point for which talairach coordinate'};
paramsInfo{end+1} = {'thisPoint',[x y z],'editable=0','The current selected point'};
for i = 1:length(gTalairach.refPoints)
  paramsInfo{end+1} = {gTalairach.refPoints{i},gTalairach.(gTalairach.refPoints{i}),sprintf('Talairach point %s',gTalairach.refPoints{i})};
end

% put up the dialog
params = mrParamsDialog(paramsInfo);

% set the point if called for
if ~isempty(params)
  % change all the rest of the points if user inputed numbers
  for i = 1:length(gTalairach.refPoints)
    gTalairach.(gTalairach.refPoints{i}) = params.(gTalairach.refPoints{i});
  end
  % set the current point to whatever user wanted
  if ~strcmp(params.whichPoint,'None')
    gTalairach.(params.whichPoint) = [x y z];
  end

  % change the center of the slice to the AC point
  gTalairach.params.xCenter = gTalairach.AC(1);
  gTalairach.params.yCenter = gTalairach.AC(2);
  gTalairach.params.zCenter = gTalairach.AC(3);

  % set the yz rotation angle
  opposite = gTalairach.PC(3) - gTalairach.AC(3);
  adjacent = gTalairach.PC(2) - gTalairach.AC(2);
  hypotenuse = sqrt(opposite^2+adjacent^2);
  gTalairach.params.yzRot = 90+r2d(acos(opposite/hypotenuse));
  
  % set the coordinates of the ACPC slice
  gTalairach.ACPCSliceCoords = calcACPCSliceCoords(gTalairach.params);

  refreshTalairachDisplay;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mouseup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mouseUpHandler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get current talairach rotation matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rotmatrix = getRotMatrix(params,offset)

global gTalairach;
dim = gTalairach.dim;

% get the transformation marix
% rotxy
c = cos(d2r(params.xyRot));
s = sin(d2r(params.xyRot));
rotxy = [c -s 0 0;s  c 0 0;0  0 1 0;0  0 0 1];

% rotyz
c = cos(d2r(params.yzRot));
s = sin(d2r(params.yzRot));
rotyz = [1  0  0 0;0  c -s 0;0  s  c 0;0  0  0 1];

% rot xz
c = cos(d2r(params.xzRot));
s = sin(d2r(params.xzRot));
rotxz = [c  0 -s 0;0  1  0 0;s  0  c 0;0  0  0  1];

% offset
if ieNotDefined('offset')
    offset = [1  0  0 params.xCenter;
              0  1  0 params.yCenter;
              0  0  1 params.zCenter;
              0  0  0    1];
end

sliceOffset = [1 0 0 -params.width/2;
	       0 1 0 -params.height/2;
	       0 0 1 0;
	       0 0 0 1];

% full rotation matrix
rotmatrix = offset*rotxy*rotyz*rotxz*sliceOffset;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function endHandler

global gTalairach;
% only run this if it is not being called by someonw else
if ~gTalairach.shutdown
    gTalairach.shutdown = 1;
    for i = 1:3
        if ishandle(gTalairach.fig(i))
            close(gTalairach.fig(i));
        end
    end
    gTalairach.shutdown = 0;
end

if gTalairach.init
  for i = 1:length(gTalairach.refPoints)
    disp(sprintf('%s: %i %i %i',gTalairach.refPoints{i},gTalairach.(gTalairach.refPoints{i})(1),gTalairach.(gTalairach.refPoints{i})(2),gTalairach.(gTalairach.refPoints{i})(3)));
  end
end

gTalairach.init = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init the interrogator handler
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initHandler(filename,vol,hdr)

global gTalairach;
disp('in init handler')
% get figure handles
gTalairach = [];
gTalairach.fig(1) = smartfig('talairach1');
gTalairach.fig(2) = smartfig('talairach2');
gTalairach.fig(3) = smartfig('talairach3');
gTalairach.roiCoords = [];
gTalairach.shutdown = 0;
gTalairach.filename = filename;
gTalairach.init = 1;
% set the callbacks appropriately
for i= 1:3
    %  set(gTalairach.fig(i),'MenuBar','none');
%    set(gTalairach.fig(i),'WindowButtonMotionFcn',sprintf('talairach(''mouseMove'',1)'));
    set(gTalairach.fig(i),'WindowButtonDownFcn',sprintf('talairach_jg(''mouseDown'',1)'));
%    set(gTalairach.fig(i),'WindowButtonUpFcn',sprintf('talairach(''mouseUp'')'));
    set(gTalairach.fig(i),'DeleteFcn',sprintf('talairach_jg(''end'',1)'));
end

% set pointer to crosshairs
set(gTalairach.fig(1),'pointer','fullcrosshair');
set(gTalairach.fig(2),'pointer','fullcrosshair');
set(gTalairach.fig(3),'pointer','fullcrosshair');

% set volume
gTalairach.vol = vol;
gTalairach.hdr = hdr;
gTalairach.dim = size(vol);
gTalairach.gamma = 0.4;

gTalairach.pos(1) = round(size(vol,1)/2);
gTalairach.pos(2) = round(size(vol,2)/2);
gTalairach.pos(3) = round(size(vol,3)/2);

% compute color map
g = gray(256);
y = g;y(:,3) = 0;
c = g;c(:,1) = 0;
m = g;m(:,2) = 0;
r = g;r(:,2:3) = 0;
myColormap = [g;y;c;r;m;m];

gTalairach.params.slices = '0';
gTalairach.params.width = 256;
gTalairach.params.height = 256;;
gTalairach.params.viewSlice = 1;
gTalairach.params.dispROI = 0;
gTalairach.params.xyRot = 0;
gTalairach.params.yzRot = 0;
gTalairach.params.xzRot = 0;
gTalairach.params.xCenter = 0;
gTalairach.params.yCenter = 0;
gTalairach.params.zCenter = 0;

% compute slice coordinates
gTalairach.params.slices = eval(gTalairach.params.slices);
gTalairach.ACPCSliceCoords = calcACPCSliceCoords(gTalairach.params);

% compute each dimension coordinates
[x y] = meshgrid(1:gTalairach.dim(2),1:gTalairach.dim(3));
gTalairach.sliceCoords{1}(1,1:length(x(:))) = gTalairach.pos(1);
gTalairach.sliceCoords{1}(2,:) = x(:);
gTalairach.sliceCoords{1}(3,:) = y(:);
gTalairach.sliceMin{1}(1:3) = [gTalairach.pos(1) min(x(:)) min(y(:))];
gTalairach.sliceMax{1}(1:3) = [gTalairach.pos(1) max(x(:)) max(y(:))];
[x y] = meshgrid(1:gTalairach.dim(1),1:gTalairach.dim(3));
gTalairach.sliceCoords{2}(1,:) = x(:);
gTalairach.sliceCoords{2}(2,:) = gTalairach.pos(2);
gTalairach.sliceCoords{2}(3,:) = y(:);
gTalairach.sliceMin{2}(1:3) = [min(x(:)) gTalairach.pos(1) min(y(:))];
gTalairach.sliceMax{2}(1:3) = [max(x(:)) gTalairach.pos(1) max(y(:))];
[x y] = meshgrid(1:gTalairach.dim(1),1:gTalairach.dim(2));
gTalairach.sliceCoords{3}(1,:) = x(:);
gTalairach.sliceCoords{3}(2,:) = y(:);
gTalairach.sliceCoords{3}(3,:) = gTalairach.pos(3);
gTalairach.sliceMin{3}(1:3) = [min(x(:)) min(y(:)) gTalairach.pos(1)];
gTalairach.sliceMax{3}(1:3) = [max(x(:)) max(y(:)) gTalairach.pos(1)];

% set up cache
gTalairach.c = mrCache('init',50);

% display three windows
figure(gTalairach.fig(1));
dispVolumeSlice(1,gTalairach.pos(1));
colormap(myColormap);

figure(gTalairach.fig(2));
dispVolumeSlice(2,gTalairach.pos(2));
colormap(myColormap);

figure(gTalairach.fig(3));
dispVolumeSlice(3,gTalairach.pos(3));
colormap(myColormap);

% set the talairach ref points
gTalairach.refPoints = {'AC','PC','PPC','AAC','LAC','RAC','SAC','IAC'}
for i = 1:length(gTalairach.refPoints)
  gTalairach.(gTalairach.refPoints{i}) = [0 0 0];
end
disp('exiting init')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% redisplay everything
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function refreshTalairachDisplay

global gTalairach;

% and redisplay
for i = 1:3
    figure(gTalairach.fig(i));
    dispVolumeSlice(i,gTalairach.pos(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate ACPC slice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ACPCSliceCoords = calcACPCSliceCoords(params)

global gTalairach;

% get the rotation matrix that transforms to the
% desired talairachd coordinates
rotmatrix = getRotMatrix(params);

% get the initial coords
% [x y] = meshgrid(-params.width/2:params.width/2,-params.height/2:params.height/2);
[x y] = meshgrid(0:(round(params.width)-1),0:(round(params.height)-1));
for i = 1:length(params.slices)
    initCoords(1,:) = x(:);
    initCoords(2,:) = y(:);
    initCoords(3,:) = params.slices(i);
    initCoords(4,:) = 1;
    % compute coordinates after transformation
    ACPCSliceCoords{i} = round(rotmatrix*initCoords);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
% display slice of volume
%%%%%%%%%%%%%%%%%%%%%%%%%
function slice = dispVolumeSlice(sliceDim,sliceNum)

global gTalairach

% save slice num
gTalairach.sliceNum(sliceDim) = sliceNum;

% get the other dimensions
otherDim = setdiff([1 2 3],sliceDim);

% clear axis so we don't keep drawing over old
cla

% get slice out of volume
switch (sliceDim)
  case {1}
    slice = squeeze(gTalairach.vol(sliceNum,:,:));
  case {2}
    slice = squeeze(gTalairach.vol(:,sliceNum,:));
  case {3}
    slice = squeeze(gTalairach.vol(:,:,sliceNum));
end

% update slice coords
gTalairach.sliceCoords{sliceDim}(sliceDim,:) = sliceNum;
gTalairach.sliceMin{sliceDim}(sliceDim) = sliceNum;
gTalairach.sliceMax{sliceDim}(sliceDim) = sliceNum;

% gamma correct
smax = max(slice(:));smin = min(slice(:));
slice = (slice-smin)./(smax-smin);
slice = 256*(slice.^gTalairach.gamma);

% now check to see if there is any overlap between this slice
% and the talairach
intersectCoords = [];displayIntersectCoords = [];
sMin = gTalairach.sliceMin{sliceDim};
sMax = gTalairach.sliceMax{sliceDim};
for i = 1:length(gTalairach.params.slices)
    % this code used to use intersect, but that is hoeplessly
    % slow. Because we are always displaying cardinal views
    % we can just do min/max checking
    ACPCSliceCoords = gTalairach.ACPCSliceCoords{i}(1:3,:)';
    thisSliceCoords = ...
        ((ACPCSliceCoords(:,1) >= sMin(1)) & ...
         (ACPCSliceCoords(:,1) <= sMax(1)) & ...
         (ACPCSliceCoords(:,2) >= sMin(2)) & ...
         (ACPCSliceCoords(:,2) <= sMax(2)) & ...
         (ACPCSliceCoords(:,3) >= sMin(3)) & ...
         (ACPCSliceCoords(:,3) <= sMax(3)));
    thisIntersectCoords = ACPCSliceCoords(thisSliceCoords,:);
    %  thisIntersectCoords = intersect(ACPCSliceCoords,gTalairach.sliceCoords{sliceDim}(1:3,:)','rows');
    if i ~= gTalairach.params.viewSlice
        intersectCoords = [intersectCoords;thisIntersectCoords];
    else
        displayIntersectCoords = thisIntersectCoords;
    end
end

% get the indexes in the current slice for those intersection coordinates
if ~isempty(intersectCoords)
    intersectIndexes = sub2ind(gTalairach.dim(otherDim),intersectCoords(:,otherDim(1)),intersectCoords(:,otherDim(2)));
else
    intersectIndexes = [];
end

if ~isempty(displayIntersectCoords)
    displayIntersectIndexes = sub2ind(gTalairach.dim(otherDim),displayIntersectCoords(:,otherDim(1)),displayIntersectCoords(:,otherDim(2)));
else
    displayIntersectIndexes = [];
end

% and set all the intersection coordinates to yellow
slice(intersectIndexes) = slice(intersectIndexes)/2+384;
% and those for the current slice to be cyan
slice(displayIntersectIndexes) = slice(displayIntersectIndexes)/2+640;

if gTalairach.params.dispROI
    % now get all the roi coords that intersect this slice
    roiCoords = gTalairach.roiCoords';
    if ~isempty(roiCoords)
        thisROICoords = ...
            ((roiCoords(:,1) >= sMin(1)) & ...
             (roiCoords(:,1) <= sMax(1)) & ...
             (roiCoords(:,2) >= sMin(2)) & ...
             (roiCoords(:,2) <= sMax(2)) & ...
             (roiCoords(:,3) >= sMin(3)) & ...
             (roiCoords(:,3) <= sMax(3)));
        thisROICoords = roiCoords(thisROICoords,:);
        if ~isempty(thisROICoords)
            % get the indexes
            roiIndexes = sub2ind(gTalairach.dim(otherDim),thisROICoords(:,otherDim(1)),thisROICoords(:,otherDim(2)));
            % and set them to red/magenta
            slice(roiIndexes) = slice(roiIndexes)+768;
        end
    end
end

% flipud and transpose
slice = flipud(slice');

% display
image(slice);
axis off;axis square

titles = {sprintf('%s\nsaggital (left-right)',gTalairach.filename),sprintf('%s\ncoronal (back-front)',gTalairach.filename),sprintf('%s\naxial (bottom-top)',gTalairach.filename)};
title(titles{sliceDim},'Interpreter','none');
dimLabels = 'xyz';
xlabel(dimLabels(otherDim(1)));
ylabel(dimLabels(otherDim(2)));
% flip the y-axis up-down
set(gca,'YTickLabel',flipud(get(gca,'YTickLabel')));
YLim = get(gca,'YLim');
set(gca,'YTick',YLim(2)-fliplr(get(gca,'YTick')));
axis on
%%%%%%%%%%%%%%%%%%%%%%%%%%
% get talairach cache val
%%%%%%%%%%%%%%%%%%%%%%%%%%
function cacheID = getCacheID(sliceNum)

global gTalairach;
p = gTalairach.params;

cacheID = sprintf('%0.1f_%0.1f_%0.1f_%0.1f_%0.1f_%0.1f_%i_%i_%0.1f',p.yzRot,p.xzRot,p.xyRot,p.xCenter,p.yCenter,p.zCenter,p.width,p.height,sliceNum);