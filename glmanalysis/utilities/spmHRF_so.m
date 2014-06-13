function hrf = spmHRF_so(TR,params)

% 2007 August 19 shani
% wrote this wrapper for spm_hrf so I can call it
% from within my glm code for mrLoadRet-4.5
% 2008 March 19 shani - added ability to do fixed effects, so have dif HRF for dif runs because dif subjects

if TR=='params' % allow user to ask for parameters
    hrf = {...
        {'description', 'spmHRF_so', 'comment describing the hdr model'},...
        {'rdelay', 6, 'response delay (spm_hrf)'},...
        {'udelay', 16, 'undershoot delay (spm_hrf)'},...
        {'rdispersion', 1, 'response dispersion (spm_hrf)'},...
        {'udispersion', 1, 'undershoot dispersion (spm_hrf)'},...
        {'incDeriv',0,'type=checkbox','include derivative of the hrf in the model?'},...
    };
    return
end

if length(params.rdelay) == 1
  p(1) = params.rdelay;
  p(2) = params.udelay;
  p(3) = params.rdispersion;
  p(4) = params.udispersion;
  p(5) = 6;
  p(6) = 0;
  p(7) = params.tmax;
  [hrf,p] = spm_hrf(TR,p);
  hrf = hrf*(1.5/TR); % adjusting from hrf (expt run at TR=1.5) to main expts (TR = 2);
else % need to have separate HRFS for separate subjects
  for iSub = 1:length(params.rdelay)
    p(1) = params.rdelay(iSub);
    p(2) = params.udelay(iSub);
    p(3) = params.rdispersion(iSub);
    p(4) = params.udispersion(iSub);
    p(5) = 6;
    p(6) = 0;
    p(7) = params.tmax;
    [hrfTmp,p] = spm_hrf(TR,p);
    hrfTmp = hrfTmp*(1.5/TR); % adjusting from hrf (expt run at TR=1.5) to main expts (TR = 2)
    hrf(:,iSub) = hrfTmp;
    clear hrfTmp;
  end
end

  