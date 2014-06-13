% first CD to the data folder

% clear MLR and mrDEFAULTS, then load them new

% run through all the analyses for the different models

% load DM, add to analysis, add as new analysis, and save new analysis with overwrite

clear analysis
>> analysis = viewGet(getMLRView,'analysis')

analysis = 

           curOverlay: 1
                    d: {[1x1 struct]}
                 date: '08-Apr-2008 17:02:47'
             function: 'eventRelatedGlm'
            groupName: 'ConcatBlur'
          guiFunction: 'eventRelatedGlmGUI'
        mergeFunction: 'defaultMergeParams'
                 name: 'glm_byDM_correctDelayS2plus'
             overlays: [1x6 struct]
               params: [1x1 struct]
    reconcileFunction: 'defaultReconcileParams'
                 type: 'glmAnal'

>> load /Users/shani/NYU/fMRI_data/GLM_output/DM/DM_correctDelayS2plus_detect_al>> bj.mat
>> analysis.d{1}.DM = DM;                   
>> viewSet(getMLRview,'newAnalysis',analysis);
>> saveAnalysis(getMLRView,analysis.name)
Warning: 'mrParamsDialog_glm_byDM_correctDelayS2plus_period_mat_already_exists'
 exceeds MATLAB's maximum name length of 63 characters and has been truncated
 to
 'mrParamsDialog_glm_byDM_correctDelayS2plus_period_mat_already_e'.
> In mrSetFigLoc at 26
  In mrParamsDialog>closeHandler at 540
  In mrParamsDialog>initFigure at 248
  In mrParamsDialog at 23
  In saveAnalysis at 48
(saveAnalysis) Overwriting old analysis
Saving /Volumes/LaCie/NYU/fMRI_data/Detection/detectionEvent/allSubj_detect/ConcatBlur/glmAnal/glm_byDM_correctDelayS2plus.mat...done
>> 
