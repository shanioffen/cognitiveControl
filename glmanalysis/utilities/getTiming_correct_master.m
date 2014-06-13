function getTiming_correct_master(subject,expt)

% syntax function getTiming_correct_master(subj,expt)
% subj: 'JG','DS','LM','RS','SO'
% expt: 'memory','detect','memDetect','vertical'
%
% It's useful to have a way to just have the starttimes, delaydurations, and whether 
% or not hte trial was correct, for each trial, rather than recalculate it from the DM
% for a given model

% check inputs and set defaults:
if nargin < 2
    disp('Must enter a subject and experiment')
    return
end

switch subject
  case 'allSubj'
    subjList = { 'JG','DS','LM','RS','SO'}; % in case allSubj
  otherwise
    subjList = {subject};
end
numSubj = length(subjList);

% Initialize some values
TR = 2; % TR is 2 sec
numTpnts = 120; % there runs were 125 time points, 
% but the first 5 are discarded, leaving 120 (=240sec)

saveDir = '/Users/shani/NYU/fMRI_data/GLM_output/DM/'; %where to save the design matrices
dataName = [saveDir 'exptTiming_' expt '_' subject];

timing.startTR = [];
timing.delayTR = [];
timing.correct = [];
runNum = 0;

for iSubj = 1:numSubj
  subj = subjList{iSubj};
  numDays = setNumDays(expt, subj);

  for nSession = 1:numDays % Go through all  days
    [dataDir expoData numRuns] = setDirectories(subj,expt,nSession);
    
    for runInd = 1:numRuns;   
      
      % keep track of total number of runs over both days and across all subj
      runNum = runNum+1;
      
      % get start times and delay durations from the expo files:
      startFile = [expoData num2str(runInd) '_startTR.xls'];
      delayFile = [expoData num2str(runInd) '_delaySec.xls'];
      corrFile = [expoData num2str(runInd) '_correct.xls'];
      startTRall = round(textread(startFile, '', 'headerlines',2)); % round to nearest TR
      delayTRall = ceil(ceil(textread(delayFile, '', 'headerlines',2))/TR); % round up to nearest sec, divide by TR, and round up to nearest TR;
      correct = textread(corrFile,'','headerlines',2);
      
      % figure out how many trials were shown in that run
      numTrials = length(correct);
      
      % limit all matrices to that number of trials (a trial may have started but not ended)
      startTR = startTRall(1:numTrials);
      delayTR = delayTRall(1:numTrials);
      
      % update startTR to account for which run we're in
      updateVal = (runNum-1) * numTpnts;
      startTR = startTR + updateVal;
      
      % stack for each run
      timing.startTR = [timing.startTR; startTR];
      timing.delayTR = [timing.delayTR; delayTR];
      timing.correct = [timing.correct; correct];
      clear startFile delayFile corrFile startTR delayTR startTRall delayTRall correct numTrials
    end % for runInd = 1:numRuns
  end % for nSession 
end % for iSubj

% save 
save(dataName,'timing')
  
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dataDir expoData numRuns]= setDirectories(subj,expt,nSession)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch expt
  case 'memory'
    switch subj
      case 'DS'
        if nSession == 1
          dataDir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/DS110904';
          expoData = '/Users/shani/NYU/Expo_files/data/memory/pilot5/denis_magnet_110904/textfiles/ds_exp';
          numRuns = 6;
        elseif nSession == 2
          dataDir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/DS112104';
          expoData = '/Users/shani/NYU/Expo_files/data/memory/pilot5/denis_magnet_112104/textfiles/ds_exp';
          numRuns = 8;
        end
      case 'JG'
        dataDir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/JG092205';
        expoData = '/Users/shani/NYU/Expo_files/data/memory/pilot5/justin_magnet_092205/textfiles/jg_exp';
        numRuns = 11;
      case 'LM'
        dataDir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/LM100305';
        expoData = '/Users/shani/NYU/Expo_files/data/memory/pilot5/leila_magnet_100305/textfiles/lm_exp';
        numRuns = 12;
      case 'RS'
        dataDir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/RS092205';
        expoData = '/Users/shani/NYU/Expo_files/data/memory/pilot5/riju_magnet_092205/textfiles/rs_exp';
        numRuns = 11; 
      case 'SO'
        dataDir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/SO051105';
        expoData = '/Users/shani/NYU/Expo_files/data/memory/pilot5/shani_magnet_051105/textfiles/so_exp';
        numRuns = 10;
    end % switch subj
    
    
  case 'memDetect'
    switch subj
      case 'DS'
        dataDir = '/Users/shani/NYU/fMRI_data/MemDetect/DS081105';
        expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot3-memDetect/denis_magnet_081105/textfiles/ds_exp';
        numRuns = 9;
      case 'JG'
        dataDir = '/Users/shani/NYU/fMRI_data/MemDetect/JG091305';
        expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot3-memDetect/justin_magnet_091305/textfiles/jg_exp';
        numRuns = 11;
      case 'LM'
        dataDir = '/Users/shani/NYU/fMRI_data/MemDetect/LM092805';
        expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot3-memDetect/leila_magnet_092805/textfiles/lm_exp';
        numRuns = 11;
      case 'RS'
        dataDir = '/Users/shani/NYU/fMRI_data/MemDetect/RS091205';
        expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot3-memDetect/riju_magnet_091205/textfiles/rs_exp';
        numRuns = 10;
      case 'SO'
        dataDir = '/Users/shani/NYU/fMRI_data/MemDetect/SO081205';
        expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot3-memDetect/shani_magnet_081205/textfiles/so_exp';    
        numRuns = 12;
    end % switch subj
    
    
  case 'detect'
    switch subj
      case 'DS'
        if nSession == 1
          dataDir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/DS050405';
          expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot1-detect/denis_magnet_050405/textfiles/ds_exp';
          numRuns = 7;
        elseif nSession == 2
          dataDir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/DS050505';
          expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot1-detect/denis_magnet_050505/textfiles/ds_exp';
          numRuns = 10;
        end
      case 'JG'
        dataDir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/JG090105';
        expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot1-detect/justin_magnet_090105/textfiles/jg_exp';
        numRuns = 10;
      case 'LM'
        dataDir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/LM083005';
        expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot1-detect/leila_magnet_083005/textfiles/lm_exp';
        numRuns = 10;
      case 'RS'
        if nSession == 1
          dataDir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/RS083005';
          expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot1-detect/riju_magnet_083005/textfiles/rs_exp';
          numRuns = 5; % there were 8 scans, but her head started hurting. 
        elseif nSession == 2
          dataDir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/RS092805';
          expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot1-detect/riju_magnet_092805/textfiles/rs_exp';
          numRuns = 5; 
        end
      case 'SO'
        if nSession == 1
          dataDir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/SO041905';
          expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot1-detect/shani_magnet_041905/textfiles/so_exp';
          numRuns = 10;
        elseif nSession == 2
          dataDir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/SO042105';
          expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot1-detect/shani_magnet_042105/textfiles/so_exp';
          numRuns = 10;
        end
    end % switch subj
    
    
  case 'vertical'
    switch subj
      case 'DS'
        dataDir = '/Users/shani/NYU/fMRI_data/Vertical/DS081505';
        expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot4-vertical/denis_magnet_081505/textfiles/ds_exp';
        numRuns = 10;
      case 'JG'
        dataDir = '/Users/shani/NYU/fMRI_data/Vertical/JG091505';
        expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot4-vertical/justin_magnet_091505/textfiles/jg_exp';
        numRuns = 11;
      case 'LM'
        dataDir = '/Users/shani/NYU/fMRI_data/Vertical/LM090605';
        expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot4-vertical/leila_magnet_090605/textfiles/lm_exp';
        numRuns = 12;
      case 'RS'
        dataDir = '/Users/shani/NYU/fMRI_data/Vertical/RS091505';
        expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot4-vertical/riju_magnet_091505/textfiles/rs_exp';
        numRuns = 11;
      case 'SO'
        dataDir = '/Users/shani/NYU/fMRI_data/Vertical/SO081805';
        expoData = '/Users/shani/NYU/Expo_files/data/attention/att_pilot4-vertical/shani_magnet_081805/textfiles/so_exp';
        numRuns = 11;
    end % switch subj
end % switch expt


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function numDays = setNumDays(expt, subj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

