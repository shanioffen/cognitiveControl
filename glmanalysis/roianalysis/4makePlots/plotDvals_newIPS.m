function plotDvals_newIPS
% function plotDvals
% this function reads in the output of the glm beta estimation,
% and creates bar graphs for the d-vals for both expts, across subj,
% by ROI
%
% 2009Sept23 - put rest in name of ROI bc restricting IPS by visLoc

restrict = 1;
HRFnum = 2;


ROIlist = {...
    'newLeftIPS_restCan25', 'newRightIPS_restCan25',...
    'newLeftIPS_restCan30', 'newRightIPS_restCan30', ...
    'newLeftIPS_restCan35', 'newRightIPS_restCan35' ...
          };
nR = 1; nC= 6;
numRoi = length(ROIlist);

subjList = {'allSubj','JG','DS','LM','RS','SO'};
numSubj = length(subjList);

exptList = {'memory','detect'};
numExpt = length(exptList);

model = 'correctTrialS2plus';

HRFlist = {'Can','Est','Vis'};
whichHRF = HRFlist{HRFnum};

figNum = figure;
plotNum = 0;

for iRoi = 1:numRoi
  ROI = ROIlist{iRoi};
  for iExpt = 1:numExpt
    expt = exptList{iExpt};
    for iSubj = 1:numSubj
      subj = subjList{iSubj};
      
      %load glm 
      glm = loadGLM(subj, expt, model, ROI, whichHRF);      
      %save dCorr for each expt,subj
      d(iSubj,iExpt) = glm.betas(2);
      dErr(iSubj,iExpt) = glm.betaError(2); %error bars based on residuals (JG)
      
      % clear the variable
      clear glm;
      
      % get ROI size per subj
      if iSubj>1
        sizeDir = '/Users/shani/NYU/fMRI_data/GLM_output/RoB_ROIanalysis/';
        sizeName = [sizeDir 'ROIsize_' subj '_' ROI];       
        load(sizeName)
        keepSize(iSubj-1,iRoi) = roiSize;
        clear roiSize;
      end
      
    end % iSubj
  end %iExpt
  
  % get average and SEM across subj for each expt
  meanD = nanmean(d(2:end,:));
  stErrD = nanstd(d(2:end,:))./sqrt(5);
  
  % calculate paired t-test and put on graph if significant
  [hValDif pValDif] = ttest(d(:,1),d(:,2));
  disp(sprintf('%s pVal dif %0.3f',ROI,pValDif))
  pValDifAll(iRoi) = pValDif;

  % also do t-test to see if each significantly different from 0
  [memZeroH memZeroP] = ttest(d(:,1));
  [detZeroH detZeroP] = ttest(d(:,2));
  disp(sprintf('%s pVal MEM %0.3f pval DET %0.3f',ROI,memZeroP,detZeroP))
  pValDetAll(iRoi) = detZeroP;
  pValMemAll(iRoi) = memZeroP;
  
  %make bar graph for average data
  plotNum = plotNum + 1;
  figure(figNum), subplot(nR,nC,plotNum);
  tempH = mybar(meanD, stErrD);
  lH(1) = tempH{1}; lH(2) = tempH{2};
  lS{1} = sprintf('memP %.2f',pValMemAll(iRoi));
  lS{2} = sprintf('detP %.2f',pValDetAll(iRoi));
  set(gca,'xticklabel','Mem|Det','fontname','arial','fontsize',8);
  title(sprintf('%s, pValDif = %.2f, size=%i',ROI,pValDifAll(iRoi),...
                round(mean(keepSize(:,iRoi))/27)));
  legend(lH,lS);
  ylim([-0.4 0.6])
  
  % save as ill file
  if iRoi == numRoi
    graphName = '/Users/shani/NYU/NYUdocuments/IPSgraphs/dPlot_MAIN_IPS_combined';
    switch whichHRF
      case 'Can'
        graphName = [graphName '_cannonHRF'];
      case 'Est'
        graphName = [graphName '_estHRF'];    
    end
    print(gcf,'-painters','-dill',graphName)
  end
  
  
  for iSubj = 1:numSubj
    subj = subjList{iSubj};
    % make a sep graph for each subj
    
    figure(figNum+1), subplot(nR*numSubj,nC,plotNum+(iSubj-1)*numRoi);
    tempH3 = mybar(d(iSubj,:),dErr(iSubj,:));
    if plotNum==1
      lH3(1) = tempH3{1};
      lS3{1} = sprintf('%s',subj);
      legend(lH3, lS3);
    end
    if iSubj>1
      roiSz = round(keepSize(iSubj-1,iRoi)/27);
    else
      roiSz = round(mean(keepSize(:,iRoi))/27);
    end
    clear tempH3
    title(sprintf('%s, size= %i',ROI,roiSz));
    set(gca,'xticklabel','Memory|Detect','fontname','arial','fontsize',8);   
    ylim([-0.6 0.8]);
  end
  

  if iRoi == numRoi
    graphName = ['/Users/shani/NYU/NYUdocuments/IPSgraphs/dPlot_bySubj_IPS_combined'];
    switch whichHRF
      case 'Can'
        graphName = [graphName '_cannonHRF'];
      case 'Est'
        graphName = [graphName '_estHRF'];    
    end
    print(gcf,'-painters','-dill',graphName)
  end
  
end %iRoi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function glm = loadGLM(subj, expt, model, ROI, whichHRF)
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


load(dataName)
if(exist('glm'))
  glm = glm;
else
  glm = [];
end
