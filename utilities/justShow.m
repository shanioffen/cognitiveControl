function justShow(view,overlayNum,scan,realX,realY,realS,roi) %,xBase,yBase,sBase)
baseCoords = viewGet(view,'mouseDownBaseCoords');
xBase = baseCoords(1); yBase = baseCoords(2); sBase = baseCoords(3);
base2tal = viewGet(view,'base2tal'); % keyboard
if(~isempty(base2tal)) % if there is a tal transform, use Tal coordinates in titles
  talCoords = round(base2tal * [xBase yBase sBase 1]');
  xTal = talCoords(1); yTal = talCoords(2); zTal = talCoords(3);
  xCoord = xTal; yCoord = yTal; sCoord = zTal;
else
  xCoord - xBase; yCoord = yBase; sCoord = sBase;
end

disp(sprintf('%i, %i, %i',xCoord,yCoord,sCoord));
  