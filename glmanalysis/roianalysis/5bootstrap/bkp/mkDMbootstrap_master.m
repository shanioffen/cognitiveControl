function sessionDesignMatrix = mkDMbootstrap_master(subj,expt)

% 03/23/2007
% This function generates design matrices with random trial
% orders (eg random delay period durations) for use in 
% running the analysis 1000 times (bootstrapping) so as to 
% calculate the range that the beta values would fall in
% if they were random (to have a measure of significance
% of the actual beta estimates for each subject)
%
% I skip some of the variations in the original DM script,
% since we have now settled on how to do the analysis:
% we do the 'noOL', that is, no overlap in the delay-period
% and the stimulus intervals; and we do the analysis for all
% the time points, rather than separating out the long and short
% delay periods.


boostrapDir = '/local/users/shani/matlab/wholeBrain/GLManalysis/ROIanalysis/5bootstrap/';
addpath(genpath(boostrapDir)); % so use the bootstrap versions of the code
% cd([bootstrapDir '/mkDM/'])

% check inputs and set defaults:
 if nargin < 2
    disp('Must enter a subject and experiment')
    return
 end

% Initialize some values
TR = 2; % TR is 2 sec
numTpnts = 120; % there runs were 125 time points, 
% but the first 5 are discarded, leaving 120 (=240sec)

% Set some variables for each subject
numDays = setNumDays(expt, subj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create design matrix
% initialize design matrix for each session:
sessionDesignMatrix = [];  

for nSession = 1:numDays % Go through all  days


    numRuns = setNumRuns(expt, subj, nSession); % use a subfunction to set some values


    for runInd = 1:numRuns

        % initialize each run's matrix to zeros
        designMatrix = zeros(numTpnts+100,3); % make it longer for overrun, will cut it down later
    
        % figure out how many trials were shown in that run
        numTrials = 20; % there woudn't have been more than that many in a single run 

	% generate a set of delay periods in random order:
	delayTR = rand(1,numTrials);
	delayTR = ceil(delayTR*7.5 + .5); % generate random delays ranging from 1 to 8 TRs uniformly
					  % to duplicate the delay sec = 1:16, but with TR = 2
	
	% use those times to create a set of startTRs
	% starting with 1
	startTR = 1;
	for trial = 1:numTrials
          
	  % set the duration of the delay for that trial:
	  delayVal = delayTR(trial);
	  
	  % Put a 1 in the s1 (1st) column at the startTR:
	  designMatrix(startTR,1) = 1;
	  
	  % Put a 1 in the s2 (3rd) column at the end of the trial (now model with 2 TRs):
	  % this will be startTR + delay - 1:
	  designMatrix(startTR+delayVal-1,3) = 1;
          designMatrix(startTR+delayVal, 3) = 1;
	  
	  % Put a 1 in the delay (2nd) column for everything in between
	  % if delay duration is only 2 TR, there's nothing in between
	  if delayVal>2
	    designMatrix(startTR+1:startTR+delayVal-2,2) = ones(delayVal-2,1);
	  end
	  
	  % update when the next trial will start: 
          % add the delay plus the 8TR ITI to the previous startTR:
	  startTR = startTR + delayVal + 8;
	  
	end % stop going through all the trials
	
	% cut the designMatrix to the right number of time points:
        designMatrix = designMatrix(1:numTpnts,:); % throw out any extra at end for unfinished trials

        % stack onto session design matrix:
        sessionDesignMatrix = [sessionDesignMatrix; designMatrix]; 

    end % for runInd = 1:numRuns
end % for nSession 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% subfunctions to make the code read more smoothly:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function numDays = setNumDays(expt, subj)
switch expt
    case 'memory'
	switch subj
	    case {'DS'}
		numDays = 2;
	    case {'JG','LM','RS','SO'}
		numDays = 1;
	    otherwise
		disp('unknown subject')
		return
	end % switch subj

    case 'memDetect'
	switch subj
	    case {'DS', 'LM', 'JG', 'RS', 'SO'}
		numDays = 1;
	    otherwise
		disp('unknown subject')
		return
	end % switch subj

    case 'detect'
	switch subj
	    case {'DS','RS','SO'}
		numDays = 2;
	    case {'JG','LM'}
		numDays = 1;
	    otherwise
		disp('unknown subject')
		return
	end % switch subj

    case 'vertical'
	switch subj
    	    case {'DS', 'LM', 'JG', 'RS', 'SO'}
		numDays = 1;
    	    otherwise
		disp('unknown subject')
		return
	end % switch subj
end % switch expt



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function numRuns = setNumRuns(expt, subj, nSession)

switch expt
 case 'memory'
  switch subj
   case 'DS'
    if nSession == 1
      numRuns = 6;
    elseif nSession == 2
      numRuns = 8;
    end
   case 'JG'
    numRuns = 11;
   case 'LM'
    numRuns = 12;
   case 'RS'
    numRuns = 11; 
   case 'SO'
    numRuns = 10;
  end % switch subj
  

 case 'memDetect'
  switch subj
   case 'DS'
    numRuns = 9;
   case 'JG'
    numRuns = 11;
   case 'LM'
    numRuns = 11;
   case 'RS'
    numRuns = 10;
   case 'SO'
    numRuns = 12;
  end % switch subj
  
  
 case 'detect'
  switch subj
   case 'DS'
    if nSession == 1
      numRuns = 7;
    elseif nSession == 2
      numRuns = 10;
    end
   case 'JG'
    numRuns = 10;
   case 'LM'
    numRuns = 10;
   case 'RS'
    if nSession == 1
      numRuns = 5; % there were 8 scans, but her head started hurting. 
    elseif nSession == 2
      numRuns = 5; 
    end
   case 'SO'
    if nSession == 1
      numRuns = 10;
    elseif nSession == 2
      numRuns = 10;
    end
  end % switch subj
  
  
 case 'vertical'
  switch subj
   case 'DS'
    numRuns = 10;
   case 'JG'
    numRuns = 11;
   case 'LM'
    numRuns = 12;
   case 'RS'
    numRuns = 11;
   case 'SO'
    numRuns = 11;
  end % switch subj
end % switch expt
