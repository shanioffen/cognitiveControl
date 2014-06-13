% from Farshad
function log2mat

X=[];
fm_logfsda;
for i=1:length(X)
    mylog = [];
    mylog.logdata = X{i};
    g1=(X{i}(:,3)==0);
    g2=(X{i}(:,3)==1);
    t = mylog.logdata(:,2);
    mylog.stimtimes_s = {t(g1), t(g2)};
    mylog.sortedby = 'suppression';
    save(['..' filesep 'FST_' num2str(i) '.mat'], 'mylog');
end;

session = [];
groups = [];
load ../../mrSession.mat

if ~strcmp(groups(4).name, 'Event-Related')
    error('Invalid groupname');
else
    groupno=4;
end
if length(X)~=length(groups(groupno).auxParams)
    error('Invalid logfile');
end

for i=1:length(groups(groupno).auxParams)
    groups(groupno).auxParams(i).stimFileName = ['FST_' num2str(i) '.mat'];
end

save ../../mrSession.mat session groups

    
