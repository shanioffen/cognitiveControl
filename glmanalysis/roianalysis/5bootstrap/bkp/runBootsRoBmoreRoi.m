% 3/26/07
% 4/20/07 - changed to run 10,000x, which requires saving and clearing intermittently
% This code runs the bootstrap
%
% 2008May11 - switched to run on RoB
%
% since you need the same DM for every ROI for a given subj/expt, this code
% goes through all the subj/expt combos, and calls another function
% that does the actual bootstrapping for all 10 ROIs.
%
% It saves the 10,000 bootstrapped values in a large matrix called bootstrap in boostrap.m
% that has five dimensions:
% numExpts x numSubj x numBetas x numROI x 10,000
% which is:
% 2 x 5 x 3 x 10 x 10,000
function bootResults = runBootsRoBmoreRoi
  
clear all

boostrapDir = '/Users/shani/NYU/matlab/wholeBrain/GLManalysis/ROIanalysis/5bootstrap/';
addpath(genpath(boostrapDir)); % so use the bootstrap versions of the code

numBoots = 5000; % how many times you want to boostrap
reRun = 0; % whether or not to remake the boostrap matrix vs just reading it in.

restrict = 1; % use data from restricted ROIs;
HRFnum = 2; %use estimated HRF
model = 'correctTrialS2plus';
hrfList = {'Can','Est','Vis'};
whichHRF = hrfList{HRFnum};

ROIlist = {...
           'leftDMsPCS','rightDMsPCS',...
           'leftDLsPCS','rightDLsPCS',...
           'leftIPCS','rightIPCS'...
           'leftAntCS','rightAntCS',...
           'leftCing','rightCing',...
           'leftPostIPS','rightPostIPS',...
           };

numRoi = length(ROIlist);

exptList = {'memory','detect'};

subjList = {'allSubj','DS','JG','LM','RS','SO'};

betaList = {'s1', 'd', 's2'};

bootdir = '/Users/shani/NYU/fMRI_data/GLM_output/RoB_ROIanalysis/bootstrap/';

total = length(subjList)*length(exptList);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run the boostrap 10,000 times to get the range of betas
for iExpt = 1:length(exptList)
  expt = exptList{iExpt};
  switch expt
    case{'detect'}
      datadir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/roiTseries/';
    case{'memory'}
      datadir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/roiTseries/';
  end
  
  for jSubj = 1:length(subjList)
    subj = subjList{jSubj};
      for iRoi = 1:numRoi
        ROI = ROIlist{iRoi};

        if(reRun)
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % load the HRF that will be used to convolve the DM
          hrfParams = loadHRFparams(subj,whichHRF,restrict);
          hrfParams.tmax = 30;
          hrf = spmHRF_so(2,hrfParams); % separate fxn;
          hipassfilter = gethipassfilter;
          subjRunList = setVariables(subj,expt);
          bootstrap=NaN*ones(length(betaList),length(ROIlist), numBoots);
          
          % get the ROI time course
          % this loads two variables: roiTseries and concatInfo
          roiDataName = [datadir subj '_' expt '_' ROI];
          if restrict
            roiDataName = [roiDataName '_restCan25.mat'];
          else
            roiDataName = [roiDataName '.mat'];
          end
          
          load(roiDataName)
          
          % get average ROI time course
          roiTC = mean(roiTseries,2);
          
          % initialize the matrix: 
          disp(sprintf('calc boots for %s %s %s',subj,expt,ROI))
          tic
            for k = 1:numBoots
              bootstrap(:,k) = bootstrapROB_master(subj,expt,roiTC,hrf,hipassfilter,subjRunList);
            end
          toc
          % save 
          saveName = [bootdir 'matFiles/bootstrap_' expt '_' subj '_' ROI];
          save(saveName,'bootstrap');
        else % if don't need to rerun
          load([bootdir 'matFiles/bootstrap_' expt '_' subj '_' ROI]);
        end % if rerun
        
        % now compare actual beta values to the bootstrapped values
        glm = loadGLM(subj,expt,model,ROI,whichHRF,restrict);
        betas = glm.betas;
        dCorr(iExpt,jSubj,iRoi) = betas(2); 
        dWrong(iExpt,jSubj,iRoi) = betas(5);
        
        % to make sure we don't get mixed up, keep these vals:
        roiName{iExpt,jSubj,iRoi} = ROI;
        subjName{iExpt,jSubj,iRoi} = subj;
        exptName{iExpt,jSubj,iRoi} = expt;
        
        numGreaterCorr(iExpt,jSubj,iRoi) = sum(bootstrap(2,:)>betas(2));
        numGreaterWrong(iExpt,jSubj,iRoi) = sum(bootstrap(2,:)>betas(5));
        clear bootstrap % before next subj/expt/roi
        
      end % going through ROIs
  end % going through subj
end % going through expt

% compare the difference between experiments
dDifCorr = squeeze(dCorr(1,:,:))-squeeze(dCorr(2,:,:));

for iSubj = 1:length(subjList)
  subj = subjList{iSubj};
  load([bootdir 'matFiles/bootstrapMoreRoi_memory_' subj]);
  bootsD = bootstrap; clear bootstrap
  load([bootdir 'matFiles/bootstrapMoreRoi_detect_' subj]);  
  bootsM = bootstrap; clear bootstrap;
  dif = squeeze(bootsM(2,:,:)) - squeeze(bootsD(2,:,:));
  
  for iRoi = 1:numRoi
    numGreaterDifCorr(iSubj,iRoi) = sum(dif(iRoi,:)>dDifCorr(iSubj,iRoi));
  end
  clear dif;
end

% turn into percents:
percentGreaterCorr = (numGreaterCorr/numBoots)*100;
percentGreaterWrong = (numGreaterWrong/numBoots)*100;
percentGreaterDifCorr = (numGreaterDifCorr/numBoots)*100;

% make it a little easier to read
total = length(exptList)*length(subjList);
pGreaterCorr = reshape(percentGreaterCorr,total,length(ROIlist));
pGreaterWrong = reshape(percentGreaterWrong,total,length(ROIlist));
ROIs = reshape(roiName,total,length(ROIlist));
SUBJs = reshape(subjName,total,length(ROIlist));
EXPTs = reshape(exptName,total,length(ROIlist));

bootResults.ROIlist = ROIlist; % save so know right order for tables
bootResults.dCorr = dCorr;
bootResults.dWrong = dWrong;
bootResults.dDifCorr = dDifCorr;
bootResults.numGreaterCorr = numGreaterCorr;
bootResults.numGreaterWrong = numGreaterWrong;
bootResults.numGreaterDifCorr = numGreaterDifCorr;
bootResults.percentGreaterCorr = percentGreaterCorr;
bootResults.percentGreaterWrong = percentGreaterWrong;
bootResults.percentGreaterDifCorr = percentGreaterDifCorr;
bootResults.pGreaterCorr = pGreaterCorr;
bootResults.pGreaterWrong = pGreaterWrong;
bootResults.ROIs = ROIs;
bootResults.SUBJs = SUBJs;
bootResults.EXPTs = EXPTs;
bootResults.roiName = roiName;
bootResults.subjName = subjName;
bootResults.exptName = exptName;

bootSaveName = [bootdir 'bootResultsMoreRoi'];
save(bootSaveName,'bootResults');

clear all;
mkTableMoreRoi % make sure to remake tables if changed anything in the boot
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function glm = loadGLM(subj, expt, model, ROI, whichHRF,restrict)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch expt
  case{'detect'}
    glmdir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/roiBetas/';
  case{'memory'}
    glmdir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/roiBetas/';
end

dataName = [glmdir 'glm_' subj '_' expt '_' ROI '_' model];
switch whichHRF
  case 'Can'
    dataName = [dataName '_cannonHRF'];
  case 'Est'
    dataName = [dataName '_estHRF'];    
end

if restrict
  dataName = [dataName '_rest'];
end

load(dataName)
if(exist('glm'))
  glm = glm;
else
  glm = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hrfParams = loadHRFparams(subj,whichHRF,restrict)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hrfDir = '/Users/shani/NYU/fMRI_data/GLM_output/HRF_est/parameters/';

% estimated HRF using detect and memory both, since need to use
% same HRF for both experiments.
glmHRFdir = '/Users/shani/NYU/fMRI_data/GLM_output/RoB_ROIanalysis/';

switch whichHRF
  case 'Can'
    hrfParams.rdelay = 6;
    hrfParams.udelay = 16;
    hrfParams.udispersion = 1;
    hrfParams.rdispersion = 1;
    if strcmp(subj,'allSubj')
      hrfParams.rdelay = 6*ones(5,1);
      hrfParams.udelay = 16*ones(5,1);
      hrfParams.udispersion = 1*ones(5,1);
      hrfParams.rdispersion = 1*ones(5,1);
    end
  
  case 'Vis'
    if strcmp(subj,'allSubj')
      % have to deal separately with each subject
      subjList = {'JG','DS','LM','RS','SO'};%switched order because looking on JG brain so put first
      for iSub = 1:length(subjList)
        sub = subjList{iSub};
        HRFfileName = [hrfDir 'SPM_HRF_params_' sub '_V1V2V3_restrict_21'];
        load(HRFfileName) % this loads the variable bestfit;
        hrfParams.rdelay(iSub,1) = bestfit.params(1);
        hrfParams.udelay(iSub,1) = bestfit.params(2);
        hrfParams.udispersion(iSub,1) = bestfit.params(3);
        hrfParams.rdispersion(iSub,1) = bestfit.params(4);
        clear bestfit;
      end
    else % just load for individual subject if not doing fixed effects
      HRFfileName = [hrfDir 'SPM_HRF_params_' subj '_V1V2V3_restrict_21'];
      load(HRFfileName) % this loads the variable bestfit
      hrfParams.rdelay = bestfit.params(1);
      hrfParams.udelay = bestfit.params(2);
      hrfParams.udispersion= bestfit.params(3);
      hrfParams.rdispersion = bestfit.params(4);
    end
  case 'Est'
    if strcmp(subj,'allSubj')
      % have to deal separately with each subject
      subjList = {'JG','DS','LM','RS','SO'};%switched order because looking on JG brain so put first
      for iSub = 1:length(subjList)
        sub = subjList{iSub};
        dataName = [glmHRFdir 'HRFest_' sub];    
        if restrict
          dataName = [dataName '_restCan25'];
        end
        load(dataName); % this loads variable estParams
        hrfParams.rdelay(iSub,1) = estParams.rdelay;
        hrfParams.udelay(iSub,1) = estParams.udelay;
        hrfParams.udispersion(iSub,1) = 1;
        hrfParams.rdispersion(iSub,1) = 1;
        clear glm;
      end
    else % just load for individual subject if not doing fixed effects
      dataName = [glmHRFdir 'HRFest_' subj];
      if restrict
        dataName = [dataName '_restCan25'];
      end
      load(dataName); % this loads variable estParams 
      hrfParams.rdelay = estParams.rdelay;
      hrfParams.udelay = estParams.udelay;
      hrfParams.udispersion = 1;
      hrfParams.rdispersion = 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hipassfilter = gethipassfilter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% since all subj/expt used same hipass filter, just take any
load /Users/shani/NYU/fMRI_data/Memory/eventRelated/roiTseries/JG_memory_rightAntCS.mat;
hipassfilter = concatInfo.hipassfilter;
clear roiTseries concatInfo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
function subjRunList = setVariables(subj,expt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
switch expt
  case{'detect'}
    % if doing fixed effects, need to keep track of how many runs each sub did, for HRF 
    % *** order of scans is JG, DS, LM, RS, SO ********    
    switch subj
      case 'JG'
        subjRunList = [10];
      case 'DS'
        subjRunList = [17];
      case 'LM'
        subjRunList = [10];
      case 'RS'
        subjRunList = [10];
      case 'SO'
        subjRunList = [20];
      case 'allSubj'
        subjRunList = [10 17 10 10 20];
    end
  case{'memory'}
    % if doing fixed effects, need to keep track of how many runs each sub did, for HRF 
    % *** order of scans is JG, DS, LM, RS, SO ********
    switch subj
      case 'JG'
        subjRunList = [11];
      case 'DS'
        subjRunList = [14];
      case 'LM'
        subjRunList = [12];
      case 'RS'
        subjRunList = [11];
      case 'SO'
        subjRunList = [10];
      case 'allSubj'
        subjRunList = [11 14 12 11 10];
    end
end
