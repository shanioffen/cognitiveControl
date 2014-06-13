mrLoadRet
analysis = viewGet(getMLRView,'analysis')
overlays = analysis.overlays
corrMinusIncorr = overlays(3);
corrMinusIncorr.name = 'dpCminusIC';  
corrMinusIncorr.type = 'dpCminusIC';
dpCorr = overlays(3); dpIncorr = overlays(6) 
corrMinusIncorr.data{1} = dpCorr.data{1} - dpIncorr.data{1};
corrMinusIncorr.range = [min(corrMinusIncorr.data{1}(:)) max(corrMinusIncorr.data{1}(:))];
corrMinusIncorr.clip =  [min(corrMinusIncorr.data{1}(:)) max(corrMinusIncorr.data{1}(:))];
v = viewSet(getMLRView,'newoverlay',corrMinusIncorr)
saveAnalysis(getMLRView,analysis.name)               
