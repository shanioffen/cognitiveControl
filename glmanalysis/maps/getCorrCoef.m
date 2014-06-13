% getCorrCoef
% 2008 May 17 Shani Offen
% adapted JG's getR2 to also
% get correlation coefficients and set as field in d, so can calculate
% pvalues and get FDR threshold for showing maps
%
function d = getCorrCoef(d)

% init some variables
correlation = [];

% check roi
slices = 1:d.dim(3);slicen = length(slices);
xvals = 1:d.dim(1);xvaln = length(xvals);
yvals = 1:d.dim(2);yvaln = length(yvals);
  
% turn off warnings to avoid divide by zero warning
warning('off','MATLAB:divideByZero');

nope = 0;
if ~nope
% display string
disppercent(-inf,'(getCorrCoef) Calculating corr and pvals');
for i = xvals
    disppercent(max((i-min(xvals))/xvaln,0.1));
  for j = yvals
    for k = slices
      timeseries = squeeze(d.data(i,j,k,d.volumes))';
      % subtract off column means
      colmeans = mean(timeseries);
      timeseries = timeseries - colmeans;
      % convert to percent signal change
      timeseries = 100*timeseries./colmeans;
      model = d.scm*squeeze(d.ehdr(i,j,k,1:d.nhdr,:));
      
      % get correlation coefficient and pval:
      [corrco pvalue] = corr(timeseries',model);
      correlation{j,k}(i) = corrco;
      pval{j,k}(i) = pvalue(1);
    end
  end
end

disppercent(inf);
end % nope  

if nope
disppercent(-inf,'(getr2) Calculating corr and pvals');
onesmatrix = ones(length(d.volumes),1);
for j = yvals
  disppercent(max((j-min(yvals))/yvaln,0.1));
  for k = slices
    % get the time series we are working on
    % this includes all the rows of one column from one slice
    % and all data points for each of these
    % thus the time series is a nxm matrix where each of the m columns
    % contains the n time points recording for that voxel
    timeseries = squeeze(d.data(:,j,k,d.volumes))';
    % subtract off column means
    colmeans = mean(timeseries,1);
    timeseries = timeseries - onesmatrix*colmeans;
    % convert to percent signal change
    timeseries = 100*timeseries./(onesmatrix*colmeans);
    temp = squeeze(d.ehdr(:,j,k,1:d.nhdr,:));
    model = (d.scm*temp');
    % get correlation coefficient and pval:
    [corrco pvalue] = corr(timeseries,model);
    correlation{j,k} = diag(corrco);
    pval{j,k} = diag(pvalue);
  end
end
disppercent(inf);
end % nope

disppercent(-inf,'(getCorrCoef) Reshaping matrices');
for i = xvals
  disppercent((i-min(xvals))/xvaln);
  for j = yvals
    for k = slices
      % now reshape correlation and pval  into a matrix
      d.correlation(i,j,k) = correlation{j,k}(i);
      d.pval(i,j,k) = pval{j,k}(i);
    end
  end
end


% display time took
disppercent(inf);

warning('on','MATLAB:divideByZero');
