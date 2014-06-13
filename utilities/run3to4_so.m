% 2007July16
% automate the conversion of my data from mrLoadRet-3 to mrLoadRet-4.5
function run3to4_so

subjList = {'DS','JG','LM','RS','SO'};
exptList = {'memory','detect','memDetect','vertical'};
total = length(subjList)*length(exptList);

h = waitbar(0,'Converting to mrLoadRet-4.5'), waitCount = 0;

% Run through all subjects and expts and days, and set variables
for iExpt = 1:length(exptList)
  for iSubj = 1:length(subjList)
    expt = exptList{iExpt};
    subj = subjList{iSubj};
    waitCount=waitCount+1; waitbar(waitCount/total,h,'Converting to mrLoadRet-4.5');
    
    if iExpt==2 & subj=='DS'
      continue % already done
    else

      % set the base directory
      switch expt
       case 'memory'
	baseDir = '~/1keystonelink/fMRI_data/Memory/eventRelated/';
       case 'detect'
	baseDir = '~/1keystonelink/fMRI_data/Detection/detectionEvent/';
       case 'memDetect'
	baseDir = '~/1keystonelink/fMRI_data/MemDetect/';
       case 'vertical'
	baseDir = '~/1keystonelink/fMRI_data/Vertical/';
      end % switch expt
      
      % set the subject directory
      nSessions = checkNumSessions(expt,subj);
      if nSessions > 1
	for iSess = 1:nSessions
	  dataDir = [baseDir subj '_' expt num2str(iSess)];
	  cd(dataDir)
	  mr3to4
	  mr3to4(4)
	end % switch subj
      else
	dataDir = [baseDir subj '_' expt];
	cd(dataDir)
	mr3to4
	mr3to4(4)
      end % going through session
    end % for skipping DS_detect
  end % for iSubj
end % for iExpt

close(h)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions

function nSessions = checkNumSessions(expt,subj)
% some subjs ran some expts on more than one day

switch expt
  
 case 'memory'
  switch subj
   case {'DS'}
    nSessions = 2;
   case {'JG','LM','RS','SO'}
    nSessions = 1;
  end
  
 case 'detect'
  switch subj
   case {'DS','RS','SO'}
    nSessions = 2;
   case {'JG','LM'}
    nSessions = 1;
  end
  
 case 'memDetect'
  switch subj
   case {'DS', 'LM', 'JG', 'RS', 'SO'}
    nSessions = 1;
  end
  
 case 'vertical'
  switch subj
   case {'DS', 'LM', 'JG', 'RS', 'SO'}
    nSessions = 1;
  end
  
end
