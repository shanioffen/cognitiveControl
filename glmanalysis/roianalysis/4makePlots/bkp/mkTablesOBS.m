tableDir = '/Users/shani/NYU/NYUdocuments/thesis/tables/Ch4/data/';
bootdir = '/Users/shani/NYU/fMRI_data/GLM_output/RoB_ROIanalysis/bootstrap/';
sizeDir = '/Users/shani/NYU/fMRI_data/GLM_output/RoB_ROIanalysis/';

ROIlist = {'leftAntCS','leftDLsPCS','leftDMsPCS','leftIPCS','leftPostIPS',...
           'rightAntCS','rightDLsPCS','rightDMsPCS','rightIPCS','rightPostIPS'};
numRoi = length(ROIlist);

exptList = {'memory','detect'}; numExpt = length(exptList);

subjList = {'JG','DS','LM','RS','SO','allSubj'}; numSubj = length(subjList);

total = length(subjList)*length(exptList);

% make tables of bootstrap results
load([bootdir 'bootResults'])
dlmwrite([tableDir 'allRoi_bootPvalCorrect.csv'],bootResults.pGreaterCorr,',');
dlmwrite([tableDir 'allRoi_bootPvalWrong.csv'],bootResults.pGreaterWrong,',');
dlmwrite([tableDir 'allRoi_dValCorr.csv'],reshape(bootResults.dCorr,total,length(ROIlist)),',');
dlmwrite([tableDir 'allRoi_dValWrong.csv'],reshape(bootResults.dWrong,total,length(ROIlist)),',');  
dlmwrite([tableDir 'allRoi_bootPvalDif.csv'],squeeze(bootResults.percentGreaterDifCorr),',');  
dlmwrite([tableDir 'alRoi_dValDif.csv'],squeeze(bootResults.dDifCorr),',');    

% also get r2 for each subj/expt/ROI and make a table
for iExpt = 1:numExpt
  expt = exptList{iExpt};
  for iSubj = 1:numSubj
    subj = subjList{iSubj};
    
    [epochr2 discard] = estimateBetas(subj,expt,1,2,0,0);
    % have to be careful because run rois in different order
    % and too hard to redo bootstrap right now
    tempr2(iExpt,iSubj,1:10) = epochr2([1 3 2 4 5 6 8 7 9 10]);
    
  end
end
r2 = reshape(tempr2,total,numRoi); % want in same shape and order at boot results;
dlmwrite([tableDir 'allRoi_R2.csv'],r2,',');

% also get ROIsizes for each subj/ROI and make a table

for iSubj = 1:numSubj-1
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
keepSize(numSubj,:) = mean(keepSize);
dlmwrite([tableDir 'allROI_roiSize.csv'],keepSize,',')
