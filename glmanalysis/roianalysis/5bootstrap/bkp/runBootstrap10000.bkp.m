% 3/26/07
% 4/20/07 - changed to run 10,000x, which requires saving and clearing intermittently
% This code runs the bootstrap
% since you need the same DM for every ROI for a given subj/expt, this code
% goes through all the subj/expt combos, and calls another function
% that does the actual bootstrapping for all 10 ROIs.
%
% It saves the 10,000 bootstrapped values in a large matrix called bootstrap in boostrap.m
% that has five dimensions:
% numExpts x numSubj x numBetas x numROI x 10,000
% which is:
% 4 x 5 x 3 x 10 x 10,000

clear all

boostrapDir = '/Users/shani/NYU/matlab/wholeBrain/GLManalysis/ROIanalysis/5bootstrap/';
addpath(genpath(boostrapDir)); % so use the bootstrap versions of the code

numBoots = 10000; % how many times you want to boostrap
reRun = 0; % whether or not to remake the boostrap matrix vs just reading it in.

ROIlist = {'V1V2V3_restrict','V1_restrict','V2_restrict', 'V3_restrict', 'V3AB_restrict', 'V4_restrict', 'V5_restrict', 'V7_restrict', 'LO1_restrict', 'LO2_restrict', 'VO1_restrict'}; 

exptList = {'memory','detect','memDetect','vertical'};

subjList = {'DS','JG','LM','RS','SO'};

betaList = {'s1', 'd', 's2'};

bootDir = '//Volumes/keystone/users/shani/fMRI_data/GLM_output/bootstrap/';

total = length(subjList)*length(exptList);
keyboard
cd(bootDir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run the boostrap 10,000 times to get the range of betas
if(reRun)
  
  % initialize the matrix: 

  
  h = waitbar(0,'calculating the bootstrap');
  counter = 0;
  for i = 1:length(exptList)
    expt = exptList{i};
  for j = 1:length(subjList)
      subj = subjList{j};

      bootstrap = NaN*ones(length(betaList),length(ROIlist), numBoots); % initialize
      for k = 1:numBoots
	bootstrap(:,:,k) = bootstrap_master(subj,expt);
%	counter = counter + 1;
%	waitbar(counter/(numBoots*total),h,['expt ' num2str(i) ' subj ' num2str(j) ' iter ' num2str(k)]);

      end
      % save and clear
      saveName = ['matFiles/bootstrap_' expt '_' subj];
      save(saveName,'bootstrap');
      clear bootstrap
    end
  end
  clear counter
  close(h)  
  
end

%% Now compare the actual betas from the real data
%% to the bootstrapped betas using random task order 
% need, for each expt/subj/ROI, to determine how many of the delay betas from the bootstrap
% are bigger than the actual measured delay beta, so we can determine if the real delay beta is 
% significantly above zero.


%%% load in the beta structure
load //Volumes/keystone/users/shani/fMRI_data/GLM_output/betas/HRFbyV1V2V3/betaStruct_HRFbyV1V2V3 % this loads the variable betaStruct, described below:

% The structure betaStruct has 200 elements (10 ROIs x 4 expt x 5 subj) and has the following fields:
% betaStruct.ROI (ROIs are {'V1V2V3_restrict', 'V1_restrict','V2_restrict', 'V3_restrict', 'V3AB_restrict', 'V4_restrict', 'V5_restrict', 
%                      'V7_restrict', 'LO1_restrict', 'LO2_restrict', 'VO1_restrict'})
% betaStruct.expt (expts are {'memory','detect','memDetect','vertical'})
% betaStruct.subj (subjects are {'DS','JG','LM','RS','SO'})
% betaStruct.betaName (to keep track of the beta order, which is {'percentDelay','betaDelay','betaFirst','betaSecond'})
% betaStruct.betaVal_wOL (the actual beta values for the above list of betas for DM with overlap)
% betaStruct.betaVal_noOL (the actual beta values for the above list of betas for DM without overlap)
%
% so the 1st 20 elements of betaStruct are all V1V2V3, the next 20 are all V1, etc. The first 5 are all memory, the next 5 are all detection, etc

% Set some variables
whichOL = 'noOL'; % default to no overlap
% Now analyze it expt by expt
% the complication is, the beta structure is organized by ROI-Expt-Subj, but we have to
% reload the boot for each Expt, which is a little silly, but it's easier for the bookkeeping

counter = 0;
for k = 1:length(ROIlist)
  for i = 1:length(exptList)
    expt = exptList{i};
    for j = 1:length(subjList)
      subj = subjList{j};
      loadName = ['bootstrap_' expt '_' subj];
      load(['//Volumes/keystone/users/shani/fMRI_data/GLM_output/bootstrap/matFiles/' loadName])

      % get beta delay for that expt/subj/ROI
      counter = counter+1;
      betaDelay(j,i,k) = betaStruct(counter).betaVal_noOL(2);
      betaFirst(j,i,k) = betaStruct(counter).betaVal_noOL(3);
      betaLast(j,i,k) = betaStruct(counter).betaVal_noOL(4);
      ROIname{j,i,k} = betaStruct(counter).ROI;
      subjName{j,i,k} = betaStruct(counter).subj;
      exptName{j,i,k} = betaStruct(counter).expt;
      
      % find out how many of those vals are bigger than the bDelay
      % and save into a matrix
      numGreater(j,i,k) = sum(bootstrap(2,k,:)>betaDelay(j,i,k));
      numFirstGreater(j,i,k) = sum(bootstrap(1,k,:)>betaFirst(j,i,k));
      clear bootstrap % reload for each expt and subj
    end
 end
end

% turn into percents:
percentGreater = (numGreater/numBoots)*100;
percentGreater = reshape(percentGreater,total,length(ROIlist));
percentFirstGreater = (numFirstGreater/numBoots)*100;
percentFirstGreater = reshape(percentFirstGreater,total,length(ROIlist));
% and just to double check
ROIs = reshape(ROIname,total,length(ROIlist));
SUBJs = reshape(subjName,total,length(ROIlist));
EXPTs = reshape(exptName,total,length(ROIlist));


save //Volumes/keystone/users/shani/fMRI_data/GLM_output/bootstrap/bootstrapResults ...
    numGreater percentGreater betaDelay ROIs SUBJs EXPTs numFirstGreater percentFirstGreater betaFirst betaLast


% also save a single 20x11 matrix with all the ROI results
dlmwrite([bootDir 'bootResults10000_all.csv'],percentGreater,','); 
dlmwrite([bootDir 'bootFirstResults.csv'], percentFirstGreater,',');

      
      
