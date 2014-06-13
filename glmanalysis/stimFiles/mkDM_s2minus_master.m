function mkDM_s2minus_master(subj,expt)

% 2008Apr6 - based on mkDM_correct_master, but in order to address
% the problems of model fit, I'm letting S2 take up 2 TRS, since that's
% what the data seem to support, and it makese sense, since a lot happens
% at the end of the trial and it may not all happen at once.
%
% I've also changed it to save out the design matrices, since maybe there's
% a problem in how the timing is being treated by Farshad's code that I'm not
% understanding.

% syntax function mkDM_correct_master(subj,expt)
% subj: 'SO','DS','RS','JG','LM'
% expt: 'memory','detect','memDetect','vertical'
%
% This program creates design matrices for each expt and subject. 
%
% The first step to estimating the beta values from the data is to create a
% design matrix from the Expo timing files. There will be one DM per scan,
% so about 10 per subject per experiment. These DMs will then be convolved
% with the HRF, de-meaned, and detrended, before being stacked into a big
% DM per experiment per subject and then regressed against the stacked time
% course from the experiments.
%
% This script is similar to the one used in the ROI based GLM analysis, but
% also takes into account whether the trials were correct or incorrect. 
% Thus there are four regressors: s1,d_correct,d_incorrect,s2.
%
% The other difference from the earlier script is that I no longer give the
% option to leave out the short delays (in the ROI based analysis we tried
% all variations with no real difference in results) and I also leave out the
% option for wOL (meaning that in this DM, the delay starts in the TR after 
% the first stimulus, and ends in the TR before the 2nd stimulus.)



% check inputs and set defaults:
if nargin < 2
    disp('Must enter a subject and experiment')
    return
end

% Initialize some values
TR = 2; % TR is 2 sec
numTpnts = 120; % there runs were 125 time points, 
% but the first 5 are discarded, leaving 120 (=240sec)

saveDir = '/Users/shani/NYU/fMRI_data/GLM_output/DM/'; %where to save the design matrices
stimDir = setStimDir(subj,expt); % where to save the stimFiles for mrLR4.5
dataName1 = [saveDir 'DM_correctDelayS2minus' expt '_' subj];
dataName2 = [saveDir 'DM_correctTrialS2minus_' expt '_' subj];

numDays = setNumDays(subj,expt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create design matrix

% initialize design matrix for each session:
DM1 = [];
DM2 = [];

for nSession = 1:numDays % Go through all  days
  [dataDir expoData numRuns] = setDirectories(subj,expt,nSession);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % first treat delays differently for correct/incorrect
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  for runInd = 1:numRuns;            
    % initialize each run's matrix to zeros
    designMatrix = zeros(numTpnts,4); % as many rows as time points (120), and four columns:
                                      % 1st stim, delay-correct, delay-incorr, 2nd stim 
    
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
    % if the subject gave the correct response, the third column
    % with a series of ones throught the delay if incorrect, 
    % and the fourth column with single 1's for the second stim
    % but no overlap between delay and S1 and S2

    for trial = 1:numTrials
      % check the delay and whether or not correct
      delayVal = delayTR(trial);
      corrVal = correct(trial);
      
      % first column: stimulus onset
      startInd = startTR(trial);
      designMatrix(startInd,1) = 1;
        
      % second and third coluns: delay
      if delayVal>3 % if delay only 1 or 2 TR, drowned out by S1 and S2
      % second column: delay when correct
        if corrVal % = 1 for correct, = 0 for incorrect
          for time = 1:delayVal-3
            designMatrix(startInd+time,2)=1;
          end 
        else % if incorrect
          for time = 1:delayVal-3
            designMatrix(startInd+time,3)=1;
          end
        end % checking if correct
      end %if delayVal>2
            
      % fourth column: second stimulus
      % secondStim happens in same TR as delay end.
      % *** Now letting it go for 3 TRs *****
      designMatrix(max(startInd,startInd+delayVal-2),4) = 1;
      designMatrix(startInd+delayVal-1,4) = 1;
      designMatrix(startInd+delayVal,4) = 1;        
    end % for trial = 1:numTrials
    
    designMatrix = designMatrix(1:120,:); % throw out any extra at end for unfinished trials
    
    % stack onto session design matrix:
    DM1 = [DM1; designMatrix];
    clear designMatrix
    clear startFile delayFile corrFile startTR delayTR correct numTrials
  end % for runInd = 1:numRuns

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Now create DM where treat the trials differently depending on whether correct or incorrect:
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  for runInd = 1:numRuns;            
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
        if delayVal>3 % if delay only 1 or 2 TR, drowned out by S1 and S2
          for time = 1:delayVal-3
            designMatrix(startInd+time,2)=1;
          end 
        end %if delayVal>2
            
        % third column: second stimulus
        % secondStim happens in same TR as delay end.
        % *** Now letting it go for 2 TRs *****
        designMatrix(max(startInd,startInd+delayVal-2),3) = 1;
        designMatrix(startInd+delayVal-1,3) = 1;
        designMatrix(startInd+delayVal,3) = 1;        
        
      else % if incorrect
        % fourth column: stimulus onset
        designMatrix(startInd,4) = 1;
        
        % fifth column: delay
        if delayVal>3 % if delay only 1 or 2 TR, drowned out by S1 and S2
          for time = 1:delayVal-3
            designMatrix(startInd+time,5)=1;
          end 
        end %if delayVal>2
            
        % sixth column: second stimulus
        % secondStim happens in same TR as delay end.
        % *** Now letting it go for 2 TRs *****   
        designMatrix(max(startInd,startInd+delayVal-2),6) = 1;
        designMatrix(startInd+delayVal-1,6) = 1;
        designMatrix(startInd+delayVal,6) = 1;        
      
      end % checking if correct
    end % for trial = 1:numTrials
    
    designMatrix = designMatrix(1:120,:); % throw out any extra at end for unfinished trials
    
    % stack onto session design matrix:
    DM2 = [DM2; designMatrix];
    clear designMatrix
    clear startFile delayFile corrFile startTR delayTR correct numTrials
  end % for runInd = 1:numRuns
end % for nSession 

% save both DMs
DM = DM1; save(dataName1, 'DM'); clear DM
DM = DM2; save(dataName2, 'DM'); clear DM

% convert to stimfiles for mrLoadRet GLM analysis:
nJunkFrames = 5;
totalRuns = length(DM1(:,1))/numTpnts;
dummyVec = (0:TR:TR*(numTpnts-1)) + (TR*nJunkFrames*(ones(1,numTpnts)));
dummyVec = dummyVec'; dummyVec = repmat(dummyVec,totalRuns,6);
timeMatrix1 = DM1.*dummyVec(:,1:4);
timeMatrix2 = DM2.*dummyVec;
for iRun = 1:totalRuns
  % get rid of extra 0's
  mylog.stimtimes_s{1}=unique(timeMatrix1(1+(iRun-1)*numTpnts:iRun*numTpnts,1)); 
  mylog.stimtimes_s{1}=mylog.stimtimes_s{1}(2:end); % if start at 0, leave this out;

  
  mylog.stimtimes_s{2}=unique(timeMatrix1(1+(iRun-1)*numTpnts:iRun*numTpnts,2)); 
  mylog.stimtimes_s{2}=mylog.stimtimes_s{2}(2:end);
  
  mylog.stimtimes_s{3}=unique(timeMatrix1(1+(iRun-1)*numTpnts:iRun*numTpnts,3)); 
  mylog.stimtimes_s{3}=mylog.stimtimes_s{3}(2:end);
  
  mylog.stimtimes_s{4}=unique(timeMatrix1(1+(iRun-1)*numTpnts:iRun*numTpnts,4)); 
  mylog.stimtimes_s{4}=mylog.stimtimes_s{4}(2:end);
  
  stimFile = [stimDir 'stimFile_corrDelay_S2minus_' num2str(iRun) '.mat'];
  save(stimFile,'mylog'), clear mylog stimFile;
  
  mylog.stimtimes_s{1}=unique(timeMatrix2(1+(iRun-1)*numTpnts:iRun*numTpnts,1)); 
  mylog.stimtimes_s{1}=mylog.stimtimes_s{1}(2:end); % if start at 0, leave this out
  
  mylog.stimtimes_s{2}=unique(timeMatrix2(1+(iRun-1)*numTpnts:iRun*numTpnts,2)); 
  mylog.stimtimes_s{2}=mylog.stimtimes_s{2}(2:end);
  
  mylog.stimtimes_s{3}=unique(timeMatrix2(1+(iRun-1)*numTpnts:iRun*numTpnts,3)); 
  mylog.stimtimes_s{3}=mylog.stimtimes_s{3}(2:end);
  
  mylog.stimtimes_s{4}=unique(timeMatrix2(1+(iRun-1)*numTpnts:iRun*numTpnts,4)); 
  mylog.stimtimes_s{4}=mylog.stimtimes_s{4}(2:end);
  
  mylog.stimtimes_s{5}=unique(timeMatrix2(1+(iRun-1)*numTpnts:iRun*numTpnts,5)); 
  mylog.stimtimes_s{5}=mylog.stimtimes_s{5}(2:end);
  
  mylog.stimtimes_s{6}=unique(timeMatrix2(1+(iRun-1)*numTpnts:iRun*numTpnts,6)); 
  mylog.stimtimes_s{6}=mylog.stimtimes_s{6}(2:end);
  
  stimFile = [stimDir 'stimFile_corrTrial_S2minus_' num2str(iRun) '.mat'];
  save(stimFile,'mylog'), clear mylog
end  % for iRun
  

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
