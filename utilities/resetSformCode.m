% had a problem when importing and concatenating, that sform_code was being reset to 1
% this code is a temporary thing to fix the sform_code by changing both the mrSession 
% and the nifti headers
% 2008Feb18 so

% code assumes you're already running a mrLR session in the folder you need to fix, since I'm
% too lazy to write the proper code to avoid doing it this way since this is just temporary
% now go through and update the headers
for iGroup = 1:viewGet(getMLRView, 'numberofGroups');
  for iScan = 1:viewGet(getMLRView, 'nScans', iGroup);
    % load the nifti header from the mrSession file
    curhdr = viewGet(getMLRView, 'niftiHdr', iScan, iGroup);
    filename = viewGet(getMLRView, 'tseriesPath', iScan, iGroup);
    % check if it is there
    if isfile(filename)
      % load the header
      hdr = cbiReadNiftiHeader(filename);
      % set the sform_code
      hdr.sform_code = 3;
      % and write it back
      hdr = cbiWriteNiftiHeader(hdr,filename);
      % read it back, (I think there is a slight numerical
      % difference in the sform44 from when it is
      % written to when it is read). This doesn't
      % affect anything, but gives better consistency
      % checking for mrUpdateNiftiHeader
      hdr = cbiReadNiftiHeader(filename);
      % now save it in the session
      v = viewSet(getMLRView,'niftiHdr',hdr,iScan,iGroup);
    else
      disp(sprintf('(resetSformCode) Could not open file %s',filename));
    end % end checking if file is there

    % and fix in the mrSession variable
    viewSet(getMLRView,'scansformcode',3,iScan,iGroup);
  end % for iScan
end % for iGroup

% save the session
saveSession

