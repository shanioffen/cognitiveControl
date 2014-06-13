function [epochR2 betasROI] = plotCorrIncorr_moreRoi(subj,expt,restrict)
% function estimateBetas(subj,expt,[restrict],[HRFnum],[redo],[doplot])
% calculate betas for given Subj,Expt for all ROIs
% plots time courses
% Model chosen is S2plus (modeling the end of hte trial with 2TRs)
%
% can use canonical HRF or estimated HRF.
% for canonical, enter HRFnum==1 or leave blank;
% for estimated, enter HRFnum==2.
% to use estimated HRF, need to first estimate it using code 'estimateHRFs.m'
  

if nargin <3
  restrict = 1;
end

ROIlist = {'leftAntCS','leftDMsPCS','leftDLsPCS','leftIPCS',...
           'rightAntCS','rightDMsPCS','rightDLsPCS','rightIPCS',...
           'leftCing','leftSupBig','leftInfBig','leftPostIPS'...
           'rightCing','rightSupBig','rightInfBig','rightPostIPS'};

numRoi = length(ROIlist);

graphDir = '/Users/shani/NYU/NYUdocuments/thesis/figures/Ch4/data/1moreRoi/';

switch expt
  case{'detect'}
    datadir = '/Users/shani/NYU/fMRI_data/Detection/detectionEvent/roiTseries/';
  case{'memory'}
    datadir = '/Users/shani/NYU/fMRI_data/Memory/eventRelated/roiTseries/';
end

% get the start TR, delayDuration, and correctStatus for calculating TTA
exptTiming = getExptTiming(subj,expt); %subfunction
figNum = figure;

plotNum = 0;
for iRoi = 1:numRoi
  ROI = ROIlist{iRoi};
  % load ROI time course and concatInfo
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
  
  % calculate TTA for plotting:
  dataTTAcorrect = getTTA(roiTC,exptTiming,1); % get TTA for correct
  dataTTAwrong = getTTA(roiTC, exptTiming,-1);
  
  plotNum = plotNum + 1; % goes across, then down
  figure(figNum), subplot(4,4,plotNum);
  legendHandle(1) = errorbar(-4:2:30, dataTTAcorrect(:,1),dataTTAcorrect(:,2),'r.-');
  legendStr{1} = 'Correct';
  hold on;
  legendHandle(2) = errorbar(-4:2:30,dataTTAwrong(:,1),dataTTAwrong(:,2),'k.-');    
  legendStr{2} = 'Incorrect';
  xlim([-4 30]);
  ylim([-0.4 0.8])
  title(sprintf('%s - %s - %s',ROI,subj,expt));
  
end %going through rois

graphName = [graphDir expt '_dataPlots_corrIncorr' subj];
print(gcf,'-painters','-dill',graphName);

