function mkDM_split_s2plus(subj,expt)


% check inputs and set defaults:
if nargin < 2
    disp('Must enter a subject and experiment')
    return
end

% Initialize some values
TR = 2; % TR is 2 sec
numTpnts = 120; % there runs were 125 time points, 
% but the first 5 are discarded, leaving 120 (=240sec)

% save 2 DMs: one where just treat the delays dif for correct/incorrect,
% the other where treat the entire trial differently for correct/incorrect:
saveDir = '/Users/shani/NYU/fMRI_data/GLM_output/DM/'; %where to save the design matrices
stimDir = setStimDir(subj,expt); % where to save the stimFiles for mrLR4.5
dataName = [saveDir 'DM_SPLITcorrectDelayS2plus_' expt '_' subj];


numDays = setNumDays(subj,expt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create design matrix

% initialize design matrix for each session:
DModd = [];
DMeven = [];

for nSession = 1:numDays % Go through all  days
  [dataDir expoData numRuns evenList oddList] = setDirectories(subj,expt,nSession);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % first DM for odd runs
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  for runInd = oddList;    % first do odd runs        
    % initialize each run's matrix to zeros
    designMatrix = zeros(numTpnts,6); % as many rows as time points (120), and six columns:
                                      % 1st stim, delay, 2nd stim for correct and incorrect
    
    % get start times and delay durations from the expo files:
    startFile = [expoData num2str(runInd) '_startTR.xls'];
    delayFile = [expoData num2str(runInd) '_delaySec.xls'];
    corrFile = [expoData num2str(runInd) '_correct.xls'];
    startTR = round(textread(startFile, '', 'headerlines',2)); % round to nearest TR
    delayTR = ceil(ceil(textread(delayFile, '', 'headerlines',2))/TR); % round up to nearest sec, divide by TR, and round up to nearest TR;
    correct = textread(corrFile,'','headerlines',2);
    
    % figure out how many trials were shown in that run
    numTrials = length(correct);
    
    % fill the first column with single 1's at every onset,
    % the second column with a series of ones throughout the delay
    % and the third column with single 1's for the second stim
    % for the correct trials; 4th, 5th and 6th columns are same for incorrect trials
    
    for trial = 1:numTrials
      % check the delay and whether or not correct
      delayVal = delayTR(trial);
      corrVal = correct(trial);
      startInd = startTR(trial);

      if corrVal % = 1 for correct
        % first column: stimulus onset
        designMatrix(startInd,1) = 1;
        
        % second column: delay
        if delayVal>2 % if delay only 1 or 2 TR, drowned out by S1 and S2
          for time = 1:delayVal-2
            designMatrix(startInd+time,2)=1;
          end 
        end %if delayVal>2
            
        % third column: second stimulus
        % secondStim happens in same TR as delay end.
        % *** Now letting it go for 2 TRs *****
        designMatrix(startInd+delayVal-1,3) = 1;
        designMatrix(startInd+delayVal,3) = 1;        
        
      else % if incorrect
        % fourth column: stimulus onset
        designMatrix(startInd,4) = 1;
        
        % fifth column: delay
        if delayVal>2 % if delay only 1 or 2 TR, drowned out by S1 and S2
          for time = 1:delayVal-2
            designMatrix(startInd+time,5)=1;
          end 
        end %if delayVal>2
            
        % sixth column: second stimulus
        % secondStim happens in same TR as delay end.
        % *** Now letting it go for 2 TRs *****        
        designMatrix(startInd+delayVal-1,6) = 1;
        designMatrix(startInd+delayVal,6) = 1;        
      
      end % checking if correct
    end % for trial = 1:numTrials
    
    designMatrix = designMatrix(1:120,:); % throw out any extra at end for unfinished trials
    
    % stack onto session design matrix:
    DModd = [DModd; designMatrix];
    clear designMatrix
    clear startFile delayFile corrFile startTR delayTR correct numTrials
  
  end % for runInd = 1:numRuns

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Now DM for even runs
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  for runInd = evenList;            
    % initialize each run's matrix to zeros
    designMatrix = zeros(numTpnts,6); % as many rows as time points (120), and six columns:
                                      % 1st stim, delay, 2nd stim for correct and incorrect
    
    % get start times and delay durations from the expo files:
    startFile = [expoData num2str(runInd) '_startTR.xls'];
    delayFile = [expoData num2str(runInd) '_delaySec.xls'];
    corrFile = [expoData num2str(runInd) '_correct.xls'];
    startTR = round(textread(startFile, '', 'headerlines',2)); % round to nearest TR
    delayTR = ceil(ceil(textread(delayFile, '', 'headerlines',2))/TR); % round up to nearest sec, divide by TR, and round up to nearest TR;
    correct = textread(corrFile,'','headerlines',2);
    
    % figure out how many trials were shown in that run
    numTrials = length(correct);
    
    % fill the first column with single 1's at every onset,
    % the second column with a series of ones throughout the delay
    % and the third column with single 1's for the second stim
    % for the correct trials; 4th, 5th and 6th columns are same for incorrect trials
    
    for trial = 1:numTrials
      % check the delay and whether or not correct
      delayVal = delayTR(trial);
      corrVal = correct(trial);
      startInd = startTR(trial);

      if corrVal % = 1 for correct
        % first column: stimulus onset
        designMatrix(startInd,1) = 1;
        
        % second column: delay
        if delayVal>2 % if delay only 1 or 2 TR, drowned out by S1 and S2
          for time = 1:delayVal-2
            designMatrix(startInd+time,2)=1;
          end 
        end %if delayVal>2
            
        % third column: second stimulus
        % secondStim happens in same TR as delay end.
        % *** Now letting it go for 2 TRs *****
        designMatrix(startInd+delayVal-1,3) = 1;
        designMatrix(startInd+delayVal,3) = 1;        
        
      else % if incorrect
        % fourth column: stimulus onset
        designMatrix(startInd,4) = 1;
        
        % fifth column: delay
        if delayVal>2 % if delay only 1 or 2 TR, drowned out by S1 and S2
          for time = 1:delayVal-2
            designMatrix(startInd+time,5)=1;
          end 
        end %if delayVal>2
            
        % sixth column: second stimulus
        % secondStim happens in same TR as delay end.
        % *** Now letting it go for 2 TRs *****        
        designMatrix(startInd+delayVal-1,6) = 1;
        designMatrix(startInd+delayVal,6) = 1;        
      
      end % checking if correct
    end % for trial = 1:numTrials
    
    designMatrix = designMatrix(1:120,:); % throw out any extra at end for unfinished trials
    
    % stack onto session design matrix:
    DMeven = [DMeven; designMatrix];
    clear designMatrix
    clear startFile delayFile corrFile startTR delayTR correct numTrials
  end % for runInd = 1:numRuns
end % for nSession 

% save both DMs
save(dataName, 'DModd','DMeven'); 

%%%%%%%%%%%%%%%%% SUBFUNCTIONS %%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function numDays = setNumDays(subj,expt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
  case 'memDetect'
    switch subj
      case {'DS', 'LM', 'JG', 'RS', 'SO'}
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dataDir expoData numRuns evenList oddList]= setDirectories(subj,expt,nSession)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch expt
  case 'memory'
    switch subj
      case 'DS'
        if nSession == 1
          dataDir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/DS110904';
          expoData = '/Users/shani/NYU/Expo_files/data/memory/pilot5/denis_magnet_110904/textfiles/ds_exp';
          numRuns = 6;
          oddList = 1:2:numRuns;
          evenList = 2:2:numRuns;
        elseif nSession == 2
          dataDir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/DS112104';
          expoData = '/Users/shani/NYU/Expo_files/data/memory/pilot5/denis_magnet_112104/textfiles/ds_exp';
          numRuns = 8;
          oddList = 1:2:numRuns;
          evenList = 2:2:numRuns;
        end
      case 'JG'
        dataDir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/JG092205';
        expoData = '/Users/shani/NYU/Expo_files/data/memory/pilot5/justin_magnet_092205/textfiles/jg_exp';
        numRuns = 11;
          oddList = 1:2:numRuns;
          evenList = 2:2:numRuns;
      case 'LM'
        dataDir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/LM100305';
        expoData = '/Users/shani/NYU/Expo_files/data/memory/pilot5/leila_magnet_100305/textfiles/lm_exp';
        numRuns = 12;
          oddList = 1:2:numRuns;
          evenList = 2:2:numRuns;
      case 'RS'
        dataDir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/RS092205';
        expoData = '/Users/shani/NYU/Expo_files/data/memory/pilot5/riju_magnet_092205/textfiles/rs_exp';
        numRuns = 11; 
          oddList = 1:2:numRuns;
          evenList = 2:2:numRuns;
      case 'SO'
        dataDir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/SO051105';
        expoData = '/Users/shani/NYU/Expo_files/data/memory/pilot5/shani_magnet_051105/textfiles/so_exp';
        numRuns = 10;
          oddList = 1:2:numRuns;
          evenList = 2:2:numRuns;
    end % switch subj
    
    
  case 'memDetect'
    switch subj
      case 'DS'
        dataDir = '/Users/shani/NYU/fMRI_data/MemDetect/DS081105';
        expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot3-memDetect/denis_magnet_081105/textfiles/ds_exp';
        numRuns = 9;
          oddList = 1:2:numRuns;
          evenList = 2:2:numRuns;
      case 'JG'
        dataDir = '/Users/shani/NYU/fMRI_data/MemDetect/JG091305';
        expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot3-memDetect/justin_magnet_091305/textfiles/jg_exp';
        numRuns = 11;
          oddList = 1:2:numRuns;
          evenList = 2:2:numRuns;
      case 'LM'
        dataDir = '/Users/shani/NYU/fMRI_data/MemDetect/LM092805';
        expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot3-memDetect/leila_magnet_092805/textfiles/lm_exp';
        numRuns = 11;
          oddList = 1:2:numRuns;
          evenList = 2:2:numRuns;
      case 'RS'
        dataDir = '/Users/shani/NYU/fMRI_data/MemDetect/RS091205';
        expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot3-memDetect/riju_magnet_091205/textfiles/rs_exp';
        numRuns = 10;
          oddList = 1:2:numRuns;
          evenList = 2:2:numRuns;
      case 'SO'
        dataDir = '/Users/shani/NYU/fMRI_data/MemDetect/SO081205';
        expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot3-memDetect/shani_magnet_081205/textfiles/so_exp';    
        numRuns = 12;
          oddList = 1:2:numRuns;
          evenList = 2:2:numRuns;
    end % switch subj
    
    
  case 'detect'
    switch subj
      case 'DS'
        if nSession == 1
          dataDir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/DS050405';
          expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot1-detect/denis_magnet_050405/textfiles/ds_exp';
          numRuns = 7;
          oddList = 1:2:numRuns;
          evenList = 2:2:numRuns;
        elseif nSession == 2
          dataDir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/DS050505';
          expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot1-detect/denis_magnet_050505/textfiles/ds_exp';
          numRuns = 10;
          oddList = 2:2:numRuns;
          evenList = 1:2:numRuns;
        end
      case 'JG'
        dataDir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/JG090105';
        expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot1-detect/justin_magnet_090105/textfiles/jg_exp';
        numRuns = 10;
          oddList = 1:2:numRuns;
          evenList = 2:2:numRuns;
      case 'LM'
        dataDir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/LM083005';
        expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot1-detect/leila_magnet_083005/textfiles/lm_exp';
        numRuns = 10;
          oddList = 1:2:numRuns;
          evenList = 2:2:numRuns;
      case 'RS'
        if nSession == 1
          dataDir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/RS083005';
          expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot1-detect/riju_magnet_083005/textfiles/rs_exp';
          numRuns = 5; % there were 8 scans, but her head started hurting. 
          oddList = 1:2:numRuns;
          evenList = 2:2:numRuns;
        elseif nSession == 2
          dataDir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/RS092805';
          expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot1-detect/riju_magnet_092805/textfiles/rs_exp';
          numRuns = 5; 
          oddList = 2:2:numRuns;
          evenList = 1:2:numRuns;
        end
      case 'SO'
        if nSession == 1
          dataDir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/SO041905';
          expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot1-detect/shani_magnet_041905/textfiles/so_exp';
          numRuns = 10;
          oddList = 1:2:numRuns;
          evenList = 2:2:numRuns;
        elseif nSession == 2
          dataDir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/SO042105';
          expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot1-detect/shani_magnet_042105/textfiles/so_exp';
          numRuns = 10;
          oddList = 1:2:numRuns;
          evenList = 2:2:numRuns;
        end
    end % switch subj
    
    
  case 'vertical'
    switch subj
      case 'DS'
        dataDir = '/Users/shani/NYU/fMRI_data/Vertical/DS081505';
        expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot4-vertical/denis_magnet_081505/textfiles/ds_exp';
        numRuns = 10;
          oddList = 1:2:numRuns;
          evenList = 2:2:numRuns;
      case 'JG'
        dataDir = '/Users/shani/NYU/fMRI_data/Vertical/JG091505';
        expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot4-vertical/justin_magnet_091505/textfiles/jg_exp';
        numRuns = 11;
          oddList = 1:2:numRuns;
          evenList = 2:2:numRuns;
      case 'LM'
        dataDir = '/Users/shani/NYU/fMRI_data/Vertical/LM090605';
        expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot4-vertical/leila_magnet_090605/textfiles/lm_exp';
        numRuns = 12;
          oddList = 1:2:numRuns;
          evenList = 2:2:numRuns;
      case 'RS'
        dataDir = '/Users/shani/NYU/fMRI_data/Vertical/RS091505';
        expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot4-vertical/riju_magnet_091505/textfiles/rs_exp';
        numRuns = 11;
          oddList = 1:2:numRuns;
          evenList = 2:2:numRuns;
      case 'SO'
        dataDir = '/Users/shani/NYU/fMRI_data/Vertical/SO081805';
        expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot4-vertical/shani_magnet_081805/textfiles/so_exp';
        numRuns = 11;
          oddList = 1:2:numRuns;
          evenList = 2:2:numRuns;
    end % switch subj
end % switch expt



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimDir = setStimDir(subj, expt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch expt
 case 'memory'
  baseDir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/';
 case 'detect'
  baseDir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/';
 case 'memDetect'
  baseDir = '/Users/shani/NYU/fMRI_data/MemDetect/';
 case 'vertical'
  baseDir = '/Users/shani/NYU/fMRI_data/Vertical/';
end

stimDir = [baseDir subj '_' expt '/Etc/stimFiles/'];

if ~isdir(stimDir)
  mkdir(stimDir)
end
