% 2007July17
% update mrSession so it knows that the stimfiles are there
% changed it to run from inside runGLM_spmHRF
% 2008March10 - changed to allow for different models
% (the difference between models:) 
% all uses 3 predictors: s1, delay, s2
% correctDelay uses 4 predictors: s1, delay_correct, delay_wrong, s2
% correctTrial uses 6 predictors: s1_correct, delay_correct, s2_correct, and s1,d,s2 _wrong
% split is like correctTrial but analyzes half the data at a time to test robustness

function stim2mrLR(expt,subj,model)

exptList = {'memory','detect','memDetect','vertical'};
subjList = {'DS','JG','LM','RS','SO','allSubj'};
total = length(subjList)*length(exptList);
h = waitbar(0,'adding stimFiles'); waitCount = 0;

waitCount=waitCount+1; waitbar(waitCount/total,h,'Adding stimFiles');

% set the directory
baseDir = setBaseDir(expt,subj);
dataDir = [baseDir subj '_' expt];
currDir = pwd;
cd(dataDir)

mrGlobals
view = newView;

groupNum = viewGet(view,'groupNum','MotionComp');
nScans = viewGet(view,'nScans',groupNum);
for iScan=2:nScans-1 % first and last scans are localizers
  stimFilename = ['stimFiles/stimFile_' model '_' num2str(iScan-1) '.mat'];
  view = viewSet(view,'stimFilename',stimFilename,iScan,groupNum);
end

saveSession;

%  end % for iSubj
%end % for iExpt

close(h);
cd(currDir)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions

%***************************************
function baseDir = setBaseDir(expt)
switch expt
 case 'memory'
  baseDir = '~/NYU/fMRI_data/Memory/eventRelated/';
 case 'detect'
  baseDir = '~/NYU/fMRI_data/Detection/detectionEvent/';
 case 'memDetect'
  baseDir = '~/NYU/fMRI_data/MemDetect/';
 case 'vertical'
  baseDir = '~/NYU/fMRI_data/Vertical/';
end % switch expt
      
