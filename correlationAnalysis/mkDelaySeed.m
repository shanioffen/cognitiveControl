% the easiest way to make the seed is to first run the GLM and then use the output
% So first open mrLR, load the GLM analysis, and then you can run this script


analysis = viewGet(getMLRView,'analysis');
delayModel = analysis.d{1}.scm(:,2);
seed{1} = 'delayActivityModel';
seed{2} = delayModel;

corrByVox(getMLRView,seed)