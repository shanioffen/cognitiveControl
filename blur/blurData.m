function blurred = blurData(filename)
% should be written as a real function
% but to start with I'm just writing it as something to be run
% from within mrLoadRet
mrGlobals;
loadMLRDefaults;

  
% read in data:
tseries = loadTSeries(getMLRView,scanNum);

% blur the data:  
blurKernel = [1 4 6 4 1]/sum([1 4 6 4 1]);
newDat = convXYZsep(data,blurKernel);

% save it back out:
% Save aveTSeries (using header of 1st scan on scanList as the template for
% the nifti header), and add it as a new scan.
scanParams.fileName = [];
scanParams.junkFrames = 0;
scanParams.nFrames = nFrames;
scanParams.description = description;
scanParams.vol2mag = vol2mag;
scanParams.vol2tal = vol2tal;
hdr = cbiReadNiftiHeader(viewGet(view,'tseriesPath',baseScan));
[viewAverage,tseriesFileName] = saveNewTSeries(viewAverage,aveTSeries,scanParams,hdr);

% Save evalstring for recomputing and params
evalstr = ['view = newView(','''','Volume','''','); view = averageTSeries(view,params);'];
[pathstr,filename,ext,versn] = fileparts(tseriesFileName);
tseriesdir = viewGet(viewAverage,'tseriesdir');
save(fullfile(tseriesdir,filename),'evalstr','params','tseriesFileName');

% Delete temporary viewBase and viewAverage
deleteView(viewBase);
deleteView(viewAverage);

return; 



  % get the path and filename
  [path,filename,ext,versn] = fileparts(viewGet(viewBase,'tseriesPath',scanNum));
  baseGroupName = viewGet(viewBase,'groupName');


    hdr = cbiReadNiftiHeader(viewGet(view,'tseriesPath',params.scanList(1)));
    % data *MUST* be written out as float32 b/c of the small values-epm
    hdr.datatype = 16;
    [viewConcat,tseriesFileName] = saveNewTSeries(viewConcat,d.data,scanParams,hdr);
    viewConcat = saveTSeries(viewConcat,d.data,saveScanNum,oldScanParams,[],1);

    % Save evalstring for recomputing and params
evalstr = ['view = newView(','''','Volume','''','); view = concatTSeries(view,params);'];
tseriesdir = viewGet(viewConcat,'tseriesdir');
[pathstr,filename,ext,versn] = fileparts(fullfile(tseriesdir,tseriesFileName));
save(fullfile(pathstr,filename),'evalstr','params','concatInfo');