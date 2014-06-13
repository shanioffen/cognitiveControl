function betas = bootstrapROB_master(subj,expt,fullRoiTC,hrf,hipassfilter,subjRunList)

% syntax: bootstrap = bootstrap_master(subj,expt)
%
% this code calls the bootstrap DM maker to get a
% set of random trial times, and then estimates the 
% betas based on those random times.
%
% it is called by runBootsRoB 10,000 per subj/expt
% runBootsRoB stores the results in a matrix with the following
% dimensions:
% beta x roi x numBoots
% 
% since you call this once for each subj/expt/ROI,
% the output of this function is a vector of size 3
% which will be fed into bigger matrix of 3 x 10 x 10000
% 
% This allows you to calculate a distribution of beta values
% for the null hypothesis, to which we can compare the 
% measured beta value.

% 04-06-07 new ROI list
% May11-2008 - changed for doing RoB, and don't get roiTcourse each time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize some values
TR = 2; % TR is 2 sec
numTpnts = 120;
betaList = {'s1', 'd', 's2'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % because some ROIs too small so no data,
  % need to get rid of NaNs for estimating betas
  % must keep track of locations so can get rid
  % of those parts of DM too 
notNanInd = find(~isnan(fullRoiTC));
if isempty(notNanInd) % means have an ROI that's empty
  roiTC = NaN*ones(length(roiTseries),1);
else
  roiTC = fullRoiTC(notNanInd);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate a random DM for the bootstrap
DMlong = mkDMboots_S2plus(subj,expt,subjRunList);
numRuns = length(DMlong)/numTpnts;

% convolve DM with the HRF to get a model matrix:
convFiltDM = [];
tempscm = [];
for scan = 1:numRuns
  % Each scan is 120 time points, so go through the matrix
  % 120 time points at a time
  DM=DMlong((scan-1)*numTpnts+1:scan*numTpnts,:);
  
  % if allsubj, determine whose runs these are
  switch subj
    case 'allSubj'
      for icount = length(subjRunList):-1:1
        if scan <= sum(subjRunList(1:icount))
          whichSub = icount;
        end
      end
    otherwise
      whichSub = 1;
  end
      
  % convolve DM with HRF
  m = convn(DM, hrf(:,whichSub));
  m = m(1:length(DM),:);
  % remove mean 
  
  m = m-repmat(mean(m), size(m,1), 1);
  % apply the same filter as original data
  if ~isempty(hipassfilter)
    m = real(ifft(fft(m) .* repmat(hipassfilter{1}', 1, size(m,2)) ));
  end
  scm = m; % no longer doing each column separately
  % stack this run's stimcmatrix on to the last one
  tempscm = [tempscm;scm];
  clear scm DM m
end

% get rid of anything corresponding to NaNs
if isempty(notNanInd) % means empty ROI
  allscm = tempscm; % keep all, will just get NaN result to save;
else
  allscm = tempscm(notNanInd,:); % save the convolved DM;
end

clear tempscm


% Now estimate the betas, quickly
precalcmatrix = ((allscm'*allscm)^-1)*allscm';
% if this don't work then do pinv
if sum(isnan(precalcmatrix(:))) == length(precalcmatrix(:))
  disp(sprintf('(estimateBetas) Using pseudo inverse to invert convolution matrix'));
  precalcmatrix = pinv(allscm);
end
    
betas = precalcmatrix*roiTC;

    

