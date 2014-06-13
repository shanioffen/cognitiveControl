function talTransform = getTalTransform(subj,points)

% given the points, produce the talairach tranformation matrix

% temporary for debugging:
% points are columns: AC, PC, PPC, AAC
% rows are dimensions: sagital, coronal, axial
switch subj
  case 'JG'

   % This is for JG
   AC = [91 140 110]';
   PC = [91 118 110]';
   PPC = [91 45 110]';
   AAC = [91 177 110]';
   LAC = [156 140 110]';
   RAC = [27 140 110]';
   SAC = [91 140 189]';
   IAC = [91 140 80]';

   % new for JG
          AC: [90 146 121]';
          PC: [90 122 113]';
         SAC: [93 122 193]';
         IAC: [87 162 75]';
         PPC: [90 45 89]';
         AAC: [90 214 144]';
         LAC: [26 145 123]';
         RAC: [154 147 119]';

   
 case 'DS'
  
  % This is for DS
  AC= [92 135 115]';
  PC= [92 108 114]';
  PPC= [92 24 110]';
  AAC= [92 206 117]';
  LAC= [29 135 114]';
  RAC= [150 135 115]';
  SAC= [92 135 186]';
  IAC= [92 135 71]';
  
 case 'LM'
    
   AC= [87 142 112]';
   PC= [87 114 103]';
   PPC= [87 46 82]';
   AAC= [87 205 131]';
   LAC= [28 142 112]';
   RAC= [151 142 112]';
   SAC= [87 142 178]';
   IAC= [93 142 69]';
   
  case 'RS'
    AC= [89 133 124]';
    PC= [89 108 115]';
    PPC= [89 44 91]';
    AAC= [89 194 145]';
    LAC= [20 133 124]';
    RAC= [156 133 124]';
    SAC= [89 133 190]';
    IAC= [89 133 81]';
    
  case 'SO'
    AC= [87 138 124]';
    PC= [87 113 116]';
    PPC= [87 39 92]';
    AAC= [87 201 143]';
    LAC= [29 138 124]';
    RAC= [149 138 124]';
    SAC= [87 138 190]';
    IAC= [87 138 81]';
    
end


points = [AC PC PPC AAC LAC RAC SAC IAC];
points(4,:) = ones(1,size(points,2));

% solve: Tal points = talTransform * original points

% Tal points are given:
tAC = [0 0 0]';
tPC = [0 -24 0]';
tPPC = [0 -102 0]';
tAAC = [0 68 0]';
tLAC = [-62 0 0]';
tRAC = [62 0 0]';
tSAC = [0 0 72]';
tIAC = [0 0 -42]';
talPoints = [tAC tPC tPPC tAAC tLAC tRAC tSAC tIAC];
talPoints(4, :) = ones(1,size(talPoints,2));

talTransform = talPoints*pinv(points);
keyboard
% can save it to a volume anatomy:
save2vol = 0; % set to 1 to save
if(save2vol)
  filename = [subj '.img'];
  newfilename = [subj 'TAL.img'];
  saveTal2baseVol(talTransform, filename, newfilename)
end

% can export it to scans
exp2scan = 0; % set to 1 to save
if(exp2scan)
  filename = [subj '.img'];
  hdr = cbiReadNiftiHeader(filename);
  exportTal2Scans(talTransform, hdr.qform44);
end


