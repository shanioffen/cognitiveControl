restrict = 1;
HRFnum = 2 
redo = 1;
doplot = 1;
ROIlist = [];
nR = [];
nC = [];


subj = 'JG';
expt = 'detect';
[epochR2 betasROI] = estimateBetas(subj,expt,restrict,HRFnum,redo,doplot,ROIlist,nR,nC);
clear epochR2 betas, close all

subj = 'JG';
expt = 'memory';
[epochR2 betasROI] = estimateBetas(subj,expt,restrict,HRFnum,redo,doplot,ROIlist,nR,nC);
clear epochR2 betas, close all;

subj = 'DS';
expt = 'detect';
[epochR2 betasROI] = estimateBetas(subj,expt,restrict,HRFnum,redo,doplot,ROIlist,nR,nC);
clear epochR2 betas, close all;

subj = 'DS';
expt = 'memory';
[epochR2 betasROI] = estimateBetas(subj,expt,restrict,HRFnum,redo,doplot,ROIlist,nR,nC);
clear epochR2 betas, close all;

subj = 'RS';
expt = 'detect';
[epochR2 betasROI] = estimateBetas(subj,expt,restrict,HRFnum,redo,doplot,ROIlist,nR,nC);
clear epochR2 betas, close all;

subj = 'RS';
expt = 'memory';
[epochR2 betasROI] = estimateBetas(subj,expt,restrict,HRFnum,redo,doplot,ROIlist,nR,nC);
clear epochR2 betas, close all;

subj = 'LM';
expt = 'detect';
[epochR2 betasROI] = estimateBetas(subj,expt,restrict,HRFnum,redo,doplot,ROIlist,nR,nC);
clear epochR2 betas, close all;

subj = 'LM';
expt = 'memory';
[epochR2 betasROI] = estimateBetas(subj,expt,restrict,HRFnum,redo,doplot,ROIlist,nR,nC);
clear epochR2 betas, close all;

subj = 'SO';
expt = 'detect';
[epochR2 betasROI] = estimateBetas(subj,expt,restrict,HRFnum,redo,doplot,ROIlist,nR,nC);
clear epochR2 betas, close all;

subj = 'SO';
expt = 'memory';
[epochR2 betasROI] = estimateBetas(subj,expt,restrict,HRFnum,redo,doplot,ROIlist,nR,nC);
clear epochR2 betas,  close all;

subj = 'allSubj';
expt = 'detect';
[epochR2 betasROI] = estimateBetas(subj,expt,restrict,HRFnum,redo,doplot,ROIlist,nR,nC);
clear epochR2 betas

subj = 'allSubj';
expt = 'memory';
[epochR2 betasROI] = estimateBetas(subj,expt,restrict,HRFnum,redo,doplot,ROIlist,nR,nC);

