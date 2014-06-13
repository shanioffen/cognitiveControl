% this function gets the FDR th if you've already loaded the analysis
% that way we can limit by more stringent THs for the maps
function getFDRth(cutoff)
analysis = viewGet(getMLRView,'analysis');
nScan = viewGet(getMLRView,'curscan');
pvals = analysis.overlays(4).data{nScan};
[th1 th2] = FDR(pvals,cutoff)