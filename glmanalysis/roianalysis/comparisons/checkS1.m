function checkS1

exptList = {'detect','memory'};
subjList = {'DS','JG','LM','RS','SO'};
roiList = {'leftDMsPCS', 'rightDMsPCS', 'leftDLsPCS', 'rightDLsPCS','newRightIPS', 'newLeftIPS'};

numExpt = length(exptList);
numSubj = length(subjList);
numRoi = length(roiList);

for iExpt = 1:numExpt
    expt = exptList{iExpt};

    for iRoi = 1:numRoi
    ROI = roiList{iRoi};

        for iSubj = 1:numSubj
        subj = subjList{iSubj};

% (1) Load data

switch expt
  case{'detect'}
    glmdir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/roiBetas/';
  case{'memory'}
    glmdir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/roiBetas/';
end

dataNameEnd = 'restFDR_correctTrialS2plus_estHRF.mat';
dataName = [glmdir 'glm_' subj '_' expt '_' ROI '_' dataNameEnd];


load(dataName)
if ~(exist('glm'))
   mrErrorDlg(['File is empty for ' expt '_' subj '_' ROI]);
   return
end

% (2) Extract the difference between S1 correct - S1 incorrect

S1dif = glm.betas(1) - glm.betas(4);



% (3) Make 2 matrices, one for Mem, one for Detect.
%     The matrix rows are subjects, columns are ROIs.

matrixName = [expt '_S1dif'];
eval([matrixName '(iSubj,iRoi) = ' num2str(S1dif) ';']);

clear glm

                end %going through subjects
       end %going through ROIs
end %going through expts


[hDet,pDet] = ttest(detect_S1dif)
[hMem,pMem] = ttest(memory_S1dif)







