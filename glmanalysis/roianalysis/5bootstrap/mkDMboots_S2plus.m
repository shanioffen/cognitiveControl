function sessionDesignMatrix = mkDMboots_S2plus(subj,expt,subjRunList)

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


% check inputs and set defaults:
 if nargin < 2
    disp('Must enter a subject and experiment')
    return
 end

% Initialize some values
TR = 2; % TR is 2 sec
numTpnts = 120; % there runs were 125 time points, but the first 5 are discarded
numRuns = sum(subjRunList);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create design matrix
% initialize design matrix for each session:
sessionDesignMatrix = [];  

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



