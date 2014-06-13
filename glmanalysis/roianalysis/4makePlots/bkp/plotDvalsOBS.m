function plotDvals(restrict,HRFnum)
% function plotDvals([restrict],[HRFnum])
% this function reads in the output of the glm beta estimation,
% and creates bar graphs for the d-vals for both expts, across subj,
% by ROI

if nargin ==0
  restrict = 1;
  HRFnum = 2;
elseif nargin == 1
  HRFnum = 2;
end


ROIlist = {'leftAntCS','leftDLsPCS','leftDMsPCS','leftIPCS','leftPostIPS',...
           'rightAntCS','rightDLsPCS','rightDMsPCS','rightIPCS','rightPostIPS'};
numRoi = length(ROIlist);

subjList = {'JG','DS','LM','RS','SO'};
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
      glm = loadGLM(subj, expt, model, ROI, whichHRF, restrict);      
      %save dCorr for each expt,subj
      d(iSubj,iExpt) = glm.betas(2);
      % clear the variable
      clear glm
      
    end % iSubj
  end %iExpt
  
  % get average and SEM across subj for each expt
  meanD = mean(d);
  stErrD = std(d)./sqrt(size(d,1));
  
  % calculate paired t-test and put on graph if significant
  [hValDif pValDif] = ttest(d(:,1),d(:,2));
  
  % also do t-test to see if each significantly different from 0
  [memZeroH memZeroP] = ttest(d(:,1));
  [detZeroH detZeroP] = ttest(d(:,2));
  
  %make bar graph
  plotNum = plotNum + 1;
  figure(figNum), subplot(2,5,plotNum);
  tempH = mybar(meanD, stErrD);
  lH(1) = tempH{1}; lH(2) = tempH{2};
  lS{1} = sprintf('%.2f',memZeroP);
  lS{2} = sprintf('%.2f',detZeroP);
  set(gca,'xticklabel','M|D','fontname','arial','fontsize',8);
  title(sprintf('%s, pValDif = %.2f',ROI,pValDif));
  legend(lH,lS);


end %iRoi

cd /Users/shani/NYU/NYUdocuments/thesis/figures/Ch4/

graphName = 'dPlot';
switch whichHRF
  case 'Can'
    graphName = [graphName '_cannonHRF'];
  case 'Est'
    graphName = [graphName '_estHRF'];    
end
if restrict
  graphName = [graphName '_rest'];
end

  
print(gcf,'-painters','-dill',graphName)


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
