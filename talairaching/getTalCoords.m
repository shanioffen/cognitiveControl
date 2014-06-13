function talVal = getTalCoords(view, varargin)

% for debugging:
if nargin<1
  talTransform = getTalTransform
  keyboard
  return
end

% given a point, return the talairach coordinates.

baseCoords = viewGet(view,'mouseDownBaseCoords');
xBase = baseCoords(1); yBase = baseCoords(2); sBase = baseCoords(3);

talTransform = getTalTransform;
orderedCoords = [xBase yBase sBase 1]';
talVal = talTransform * orderedCoords;
talVal = talVal(1:3)

function talTransform = getTalTransform
% eventually should read from s-form.

% given the points, produce the talairach tranformation matrix
% temporarily get by hand:
if nargin < 1
  % this is from doing it in mrLR
  % AC = [91 140 110]';
  % PC = [91 118 110]';
  % PPC = [91 45 110]';
  % AAC = [91 177 110]';
  % LAC = [156 140 110]';
  % RAC = [27 140 110]';
  % SAC = [91 140 189]';
  % IAC = [91 140 80]';
  
  % this is from choosing points using JG's talairach (adapted from reslice)
  AC = [90 148 121]';
  PC = [90 121 114]';
  PPC = [90 47 94]';
  AAC = [90 214 138]';
  LAC = [32 148 126]';
  RAC = [154 148 116]';
  SAC = [97 148 186]';
  IAC = [88 148 77]';
  
  % points are columns: AC, PC, PPC, AAC
  % rows are dimensions: sagital, coronal, axial
  points = [AC PC PPC AAC LAC RAC SAC IAC];
  points(4,:) = ones(1,size(points,2));
  save talPoints points
end
% solve: Tal points = talTransform * original points

% Tal points are given:
tAC = [0 0 0]';
tPC = [0 -24 0]';
tPPC = [0 -102 0]';
tAAC = [0 68 0]';
tLAC = [-62 0 0]';
tRAC = [62 0 0]';
tSAC = [0 0 72]';
tIAC = [0 0 -42]';
talPoints = [tAC tPC tPPC tAAC tLAC tRAC tSAC tIAC];
talPoints(4, :) = ones(1,size(talPoints,2));

talTransform = talPoints*pinv(points);

% save the transform
save talTransform talTransform

