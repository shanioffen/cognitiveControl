function mkDM_allSubj

% This just stacks the DMs for all the subjects in the following order:
% JG, DS, LM, RS, SO

subjList = {'JG','DS','LM','RS','SO'};
exptList = {'detect','memory'};

modelList = {'correctTrial','correctTrialS2plus','s2plusLong','correctTrialS2trip','S2tripLong'};
DMdir = '/Users/shani/NYU/fMRI_data/GLM_output/DM/'; 

for iExpt = 1:length(exptList)
  expt = exptList{iExpt};
  bigDM = [];
  
  for iModel = length(modelList)
    model = modelList{iModel};
    bigDM = [];
  
    for iSubj = 1:length(subjList)
      subj = subjList{iSubj};
    
      load([DMdir 'DM_' model '_' expt '_' subj]); % loads variable DM;
      bigDM = [bigDM; DM];
      clear DM;
    end % for iSubj
    
    DM = bigDM; % want the variable saved to be called DM;
    save([DMdir 'DM_' model '_' expt '_allSubj'],'DM');
    
  end % going through models
end % going through expts


%modelList = {'correctDelay','correctTrial','correctDelayS2plus','correctTrialS2plus','correctDelayS2minus','correctTrialS2minus'};