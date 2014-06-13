tableDir = '/Users/shani/NYU/NYUdocuments/thesis/tables/Ch4/data/1moreRoi/';
bootdir = '/Users/shani/NYU/fMRI_data/GLM_output/RoB_ROIanalysis/bootstrap/';
sizeDir = '/Users/shani/NYU/fMRI_data/GLM_output/RoB_ROIanalysis/';

ROIlist = {'leftAntCS','leftDMsPCS','leftDLsPCS','leftIPCS',...
           'rightAntCS','rightDMsPCS','rightDLsPCS','rightIPCS',...
           'leftCing','leftSupBig','leftInfBig','leftPostIPS'...
           'rightCing','rightSupBig','rightInfBig','rightPostIPS'};
clear ROIlist; % now get it from bootResults so don't confuse the order

exptList = {'memory','detect'}; numExpt = length(exptList);

subjList = {'allSubj','DS','JG','LM','RS','SO'}; numSubj = length(subjList);

total = length(subjList)*length(exptList);

% make tables of bootstrap results
load([bootdir 'bootResultsMoreRoi']);
ROIlist = bootResults.ROIlist; % so know what order were run in for boots;
numRoi = length(ROIlist);
dlmwrite([tableDir 'allRoi_bootPvalCorrect.csv'],bootResults.pGreaterCorr,',');
dlmwrite([tableDir 'allRoi_bootPvalWrong.csv'],bootResults.pGreaterWrong,',');
dlmwrite([tableDir 'allRoi_dValCorr.csv'],reshape(bootResults.dCorr,total,length(ROIlist)),',');
dlmwrite([tableDir 'allRoi_dValWrong.csv'],reshape(bootResults.dWrong,total,length(ROIlist)),',');  
dlmwrite([tableDir 'allRoi_bootPvalDif.csv'],squeeze(bootResults.percentGreaterDifCorr),',');  
dlmwrite([tableDir 'allRoi_dValDif.csv'],squeeze(bootResults.dDifCorr),',');    

% also get r2 for each subj/expt/ROI and make a table
for iExpt = 1:numExpt
  expt = exptList{iExpt};
  for iSubj = 1:numSubj
    subj = subjList{iSubj};
    
    [epochr2 discard] = estimateBetasMorerois(subj,expt,1,2,0,0,ROIlist);
    % have to be careful because run rois in different order
    % and too hard to redo bootstrap right now
    tempr2(iExpt,iSubj,:) = epochr2;
    
  end
end
r2 = reshape(tempr2,total,numRoi); % want in same shape and order at boot results;
dlmwrite([tableDir 'allRoi_R2_moreRoi.csv'],r2,',');

% also get ROIsizes for each subj/ROI and make a table
keepsize = NaN*ones(numSubj,numRoi);

for iSubj = 2:numSubj
  subj = subjList{iSubj};
  for iRoi = 1:numRoi
    ROI = ROIlist{iRoi};
    
    sizeName = [sizeDir 'ROIsize_' subj '_' ROI];
    sizeName = [sizeName '_restCan25.mat'];
    load(sizeName)
    keepSize(iSubj,iRoi) = roiSize;
    clear roiSize;
    
  end
end
keepSize(1,:) = nanmean(keepSize);
dlmwrite([tableDir 'allROI_roiSize.csv'],keepSize,',')
