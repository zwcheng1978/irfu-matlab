function [timeIntervalsOfInterest]=sc_separation_sorting(varargin)

%% find time intervals when Cluster is in tailbox
Units=irf_units;
clusterPositionFileGSE = '/data/caalocal/mR_1min.mat';
irf_log('fcal',['Cluster position file: ' clusterPositionFileGSE]);
if ~exist(clusterPositionFileGSE,'file')
    irf.log('fcal','Position file does not exist, log in to spis or brain!')
    return;
end
disp('Loading Cluster 1 min positions');
load(clusterPositionFileGSE,'R'); % load Cluster positions 1min resolution
%Only interested in data from 2003
tint=[irf_time([2003 07 01 01 00 00]) irf_time([2004 01 01 01 00 00])];
R=irf_tlim(R,tint(1),tint(2),0);

tStep=median(diff(R(:,1))); % time step
tailBoxX=-4; % tailbox is at X less than this value, in REdd
tailBoxDZ=10; % tailbox distance halfwidth in Z, in RE
tailBoxS=2;% |Y| < |X|^tailBoxS, for X=-5, |Y|<25 and for X=-20, |Y|<400
R=irf_gse2gsm(R);
RRE=R;
RRE(:,2:end)=R(:,2:end)*Units.km/Units.RE;clear R; %Calculates the s/c positions in RE units

% tailbox definition
conditionString = ['X<' num2str(tailBoxX) 'RE,|Z|<' num2str(tailBoxDZ) 'RE,|Y|<|X|^' num2str(tailBoxS) ' RE GSM'];
disp(['Finding when Cluster is in tailbox, ' conditionString]);
% all indexes when Cluster in tailbox
itailbox = (abs(RRE(:,3))<abs(RRE(:,2)).^tailBoxS & abs(RRE(:,4))<tailBoxDZ & RRE(:,2)<tailBoxX);
% put the interval start/end time half a step before/after start/end index times
indstart                = find(diff([0 itailbox(:)']) == 1)';
indend                  = find(diff([itailbox(:)' 0]) == -1)';
timeIntervalsOfInterest = [RRE(indstart,1)-tStep/2 RRE(indend,1)+tStep/2];

scSeparationFulfilled=false(length(timeIntervalsOfInterest(:,1)),1);
disp('Calculating and sorting based on s/c separation');
for iTimeInterval=1:length(timeIntervalsOfInterest(:,1))
%Calculates the separation between the satellites
tint=timeIntervalsOfInterest(iTimeInterval,:);
R1 = local.c_read('R1',tint);
R2 = local.c_read('R2',tint);
R3 = local.c_read('R3',tint);
R4 = local.c_read('R4',tint);
R1=irf_gse2gsm(R1);
R2=irf_gse2gsm(R2);
R3=irf_gse2gsm(R3);
R4=irf_gse2gsm(R4);
if isempty(R1)
    tint
    continue
end
if isempty(R2)
    tint
    continue
end
if isempty(R3)
    tint
    continue
end
if isempty(R4)
    tint
    continue
end
R2=irf_resamp(R2,R1);
R3=irf_resamp(R3,R1);
R4=irf_resamp(R4,R1);
dR1R2=irf_add(1,R1,-1,R2);
dR1R2=irf_abs(dR1R2);
dR1R3=irf_add(1,R1,-1,R3);
dR1R3=irf_abs(dR1R3);
dR1R4=irf_add(1,R1,-1,R4);
dR1R4=irf_abs(dR1R4);
dR2R3=irf_add(1,R2,-1,R3);
dR2R3=irf_abs(dR2R3);
dR2R4=irf_add(1,R2,-1,R4);
dR2R4=irf_abs(dR2R4);
dR3R4=irf_add(1,R3,-1,R4);
dR3R4=irf_abs(dR3R4);
d_i=1000;
sortDistance=max([dR1R2(:,5) dR1R3(:,5) dR1R4(:,5) dR2R4(:,5) dR2R3(:,5) dR3R4(:,5)],[],2) <= d_i;

if all(~sortDistance)   %If none of the sort is fullfilled then just continue
    continue;
end
    
% Checks if the amount of times during the interval that didn't fulfill the
% requirement is less than 1/3 of the time interval the the scseparation is
% assumed to be acceptable
Rtest=R1;
Rtest(~sortDistance,2:4)=NaN;
if sum(isnan(Rtest(:,2)))<=((1/3)*length(R1(:,2)))
scSeparationFulfilled(iTimeInterval,1)=true;
end
end

%Removes time intervals (tint) where the separation between the s/c's weren't
%acceptable
timeIntervalsOfInterest(~scSeparationFulfilled,:)=NaN;
end