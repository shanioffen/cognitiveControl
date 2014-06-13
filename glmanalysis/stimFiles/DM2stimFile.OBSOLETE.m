% 2007July16
%
% code to convert my DM into stimFiles for the GLM analysis in the new mrLoadRet4.5
% The stimFiles have to have variable name mylog, with field name stimtimes_s, 
% with the stimulus times for each condition, 
% 2007July24 - fixed a bug
%
% To do - want to make use of the stimdurations_s field
% but will need to read in the excel files, so for the meantime
% I'm just treating the delay duration as a series of consecutive events


function DM2stimFile

% List all subj and expt

subjList = {'DS','JG','LM','RS','SO'};
exptList = {'memory','detect','memDetect','vertical'};
total = length(subjList)*length(exptList);

TR = 2;
nTpnts = 120;
nJunkFrames = 5; % throw out first 5 frames

h = waitbar(0,'Converting Design Matrices'); waitCount = 0;
% Run through all subjects and expts and days, and set variables
for iExpt = 1:length(exptList)
  for iSubj = 1:length(subjList)
    expt = exptList{iExpt};
    subj = subjList{iSubj};
    
    DMname = ['DM_' expt '_' subj '_noOL_0.mat'];
    DMfile = ['~/1keystonelink/matlab/noGUI/GLM/DM/' DMname];
    stimDir = setStimDir(subj,expt); % neater to make this a subfunction

    waitCount = waitCount+1; waitbar(waitCount/total,h,['Converting DM for ' subj ' ' expt]);
    
    % Load the DM already created; this load sessionDesignMatrix, with 6 columns; we only care about the last 3 columns0
    load(DMfile)

    % calculate how many runs that subj did for that expt
    numRuns = length(sessionDesignMatrix(:,1))/nTpnts;

    % need to convert from 1's and 0's to times, by multiplying:
    % starting at first TR (2), might need to start at 0
    dummyVec = (TR:TR:TR*(nTpnts-0)) + (TR*nJunkFrames*(ones(1,nTpnts)));
    dummyVec = dummyVec'; dummyVec = repmat(dummyVec,numRuns,3);
    timeMatrix = sessionDesignMatrix(:,4:6).*dummyVec; clear dummyVec
    
    % Convert it to the stimFile format
    % run through it in 120-element blocks, since each run was 120 time points, and need a separate stimFile for each run
    
    for iRun = 1:numRuns
      mylog.stimtimes_s{1}=unique(timeMatrix(1+(iRun-1)*nTpnts:iRun*nTpnts,1)); % get rid of extra 0's
      mylog.stimtimes_s{1}=mylog.stimtimes_s{1}(2:end); % if start at 0, leave this out

      mylog.stimtimes_s{2}=unique(timeMatrix(1+(iRun-1)*nTpnts:iRun*nTpnts,2)); 
      mylog.stimtimes_s{2}=mylog.stimtimes_s{2}(2:end);
      
      mylog.stimtimes_s{3}=unique(timeMatrix(1+(iRun-1)*nTpnts:iRun*nTpnts,3)); 
      mylog.stimtimes_s{3}=mylog.stimtimes_s{3}(2:end);
      
      stimFile = [stimDir 'stimFile_' num2str(iRun) '.mat'];
      save(stimFile,'mylog'), clear mylog
    end  % for iRun
    clear timeMatrix sessionDesignMatrix

  end % for iSubj
end % for iExpt

close(h)

%%%% subfunction to set the stimulus directory for output
function stimDir = setStimDir(subj, expt)

switch expt
 case 'memory'
  baseDir = '~/1keystonelink/fMRI_data/Memory/eventRelated/';
 case 'detect'
  baseDir = '~/1keystonelink/fMRI_data/Detection/detectionEvent/';
 case 'memDetect'
  baseDir = '~/1keystonelink/fMRI_data/MemDetect/';
 case 'vertical'
  baseDir = '~/1keystonelink/fMRI_data/Vertical/';
end

stimDir = [baseDir subj '_' expt '/Etc/stimFiles/'];

if ~isdir(stimDir)
  mkdir(stimDir)
end




