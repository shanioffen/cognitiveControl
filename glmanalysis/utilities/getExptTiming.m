%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function timing = getExptTiming(subj, expt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

saveDir = '/Users/shani/NYU/fMRI_data/GLM_output/DM/'; %where to save the design matrices
dataName = [saveDir 'exptTiming_' expt '_' subj];
load(dataName);
if(exist('timing'))
  timing = timing;
else
  timing = [];
  disp(sprintf('(getExptTiming): No timing file exists for subj %s expt %s',subj, expt))
end
