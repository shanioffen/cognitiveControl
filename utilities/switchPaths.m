% mrPaths.m
%
% edited 2007Aug30 for Katrin

%      usage: mrPaths.m()
%         by: justin gardner
%       date: 03/02/05
%    purpose: load paths for mrVISTA software
%
function retval = switchPaths(ver)

if (nargin == 0)
  ver = 3.1;
end

shaniPath = '~/1keystonelink/matlab';
TFIpath = '/share/wotan/heegerlab';
heegerpath = '/share/wotan/heegerlab';

disp(sprintf('Adding mrVISTA paths for version %0.2f',ver));

% remove all possible old paths
removepath(fullfile(heegerpath,'mrUtilities'));
removepath(fullfile(heegerpath,'mrLoadRet-3.0'));
removepath(fullfile(heegerpath,'mrAlign-4.1'));
removepath(fullfile(heegerpath,'mrUtilities-3.1'));
removepath(fullfile(heegerpath,'mrLoadRet-3.1'));
removepath(fullfile(heegerpath,'mrAlign-4.2'));
removepath(fullfile(shaniPath,'mrUtilities'));
removepath(fullfile(shaniPath,'mrLoadRet-3.0'));
removepath(fullfile(shaniPath,'mrAlign-4.1'));
removepath(fullfile(shaniPath,'mrUtilities-3.1'));
removepath(fullfile(shaniPath,'mrLoadRet-3.1'));
removepath(fullfile(shaniPath,'mrAlign-4.2'));
removepath(fullfile(shaniPath,'mrTools-4.5'));

if (ver == 3.1)
  addpath(fullfile(TFIpath,'TFI/macosx/matlab'));
  addpath(genpath(fullfile(heegerpath,'mrUtilities-3.1')));
  addpath(genpath(fullfile(heegerpath,'mrLoadRet-3.1')));
elseif (ver >= 4.0)
  addpath(genpath(fullfile(shanipath,'mrTools-4.5')));
else
  disp(sprintf('UHOH: Unknown version %f',ver));
  return
end

setpref('VISTA','site','NYU');

% ****************************************** %
% CHANGE THIS LINE TO HAVE YOUR ANATOMY PATH %
% setpref('VISTA', 'defaultAnatomyPath', '/Local/Users/justin/data/nyu/');
% ****************************************** %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% only remove pathname if it exists
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function removepath(pathname)

if ~isempty(findstr(path,pathname))
  rmpath(genpath(pathname));
end
