function [timeIntervalWithinInterval, timeIntervalsOfInterest]=c_4_null_sorting
%C_4_NULL_SORTING - Sorts data for finding nulls in the magnetotail in 2003
%from July to 1 Jan 2004.
%
%This function sorts data by only including data measured in the
%magnetotail that fulfills the restrictions of maximum spacecraft
%separation of one ion inertial length and that the position of the null is
%within the spacecraft the spacecraftbox (maximum distance of one ion inertial length).
%
%   [timeIntervalWithinInterval, timeIntervalsOfInterest]=c_4_NULL_SORTING
%   
%   OUTPUT
%   timeIntervalWithinInterval: the time intervals within the time intervals where the satellite was in ...
%the magneotail (defined here as x<-4RE, y<16, z<10) where nulls were found.
%   timeIntervalsOfInterest The time intervals where 2/3 of the times had
%   nulls in it.
%   INPUT
%   Just function name
%
%See Also C_4_NULL_POSITION
%

%Removes time intervals (tint) where the separation between the s/c's weren't
%acceptable in sc_separation_sorting function

%Sorting data so that the maximum distance to a null is smaller or equal to the
%largest distance between all satellites thus keeping it to an ion
%inertial length where Taylor expansion is assumed to work (Bfield behaves
%linearly)
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

disp('Starting calculation based on the magnetic field and position of the null');
nullPositionFulfilled = false(length(timeIntervalsOfInterest(:,1)),1);
timeIntervalWithinInterval = zeros(length(timeIntervalsOfInterest(:,1)),2);
datagap=zeros(1,2);
for iTimeInterval2=1:length(timeIntervalsOfInterest(:,1))
tint=timeIntervalsOfInterest(iTimeInterval2,:); %Goes through each tint every loop
timestepstr=num2str(iTimeInterval2);
disp(['Timestep: ' timestepstr]);
if all(isnan(tint)) %If the tint has been NaN (the separation between the 
    %satellites wasn't fulfilled) continue to next tint
    continue
end
if iTimeInterval2==35
B1=local.c_read('B_vec_xyz_gse__C1_CP_FGM_5VPS',tint);
B1 = irf_gse2gsm(B1); %Gives the magnetic field of Cluster 1 in GSM
B2=local.c_read('B_vec_xyz_gse__C2_CP_FGM_5VPS',tint);
B2 = irf_gse2gsm(B2); %Gives the magnetic field of Cluster 1 in GSM
B3=local.c_read('B_vec_xyz_gse__C3_CP_FGM_5VPS',tint);
B3 = irf_gse2gsm(B3); %Gives the magnetic field of Cluster 1 in GSM
B4=local.c_read('B_vec_xyz_gse__C4_CP_FGM_5VPS',tint);
B4 = irf_gse2gsm(B4); %Gives the magnetic field of Cluster 1 in GSM    
else
B1=local.c_read('B_vec_xyz_gse__C1_CP_FGM_FULL',tint);
B1 = irf_gse2gsm(B1); %Gives the magnetic field of Cluster 1 in GSM
B2=local.c_read('B_vec_xyz_gse__C2_CP_FGM_FULL',tint);
B2 = irf_gse2gsm(B2); %Gives the magnetic field of Cluster 1 in GSM
B3=local.c_read('B_vec_xyz_gse__C3_CP_FGM_FULL',tint);
B3 = irf_gse2gsm(B3); %Gives the magnetic field of Cluster 1 in GSM
B4=local.c_read('B_vec_xyz_gse__C4_CP_FGM_FULL',tint);
B4 = irf_gse2gsm(B4); %Gives the magnetic field of Cluster 1 in GSM
end
R1 = local.c_read('R1',tint);
R2 = local.c_read('R2',tint);
R3 = local.c_read('R3',tint);
R4 = local.c_read('R4',tint);

if isempty(B1)
    datagap(iTimeInterval2,:)=tint;
    continue
end
if isempty(B2)
    datagap(iTimeInterval2,:)=tint;
    continue
end
if isempty(B3)
    datagap(iTimeInterval2,:)=tint;
    continue
end
if isempty(B4)
    datagap(iTimeInterval2,:)=tint;
    continue
end
if isempty(R1)
   datagap(iTimeInterval2)=tint;
    continue
end
if isempty(R2)
    datagap(iTimeInterval2,:)=tint;
    continue
end
if isempty(R3)
    datagap(iTimeInterval2,:)=tint;
    continue
end
if isempty(R4)
    datagap(iTimeInterval2,:)=tint;
    continue
end

R1=irf_gse2gsm(R1);
R2=irf_gse2gsm(R2);
R3=irf_gse2gsm(R3);
R4=irf_gse2gsm(R4);

%Resample B- and R-fields
B2=irf_resamp(B2,B1);
B3=irf_resamp(B3,B1);
B4=irf_resamp(B4,B1);
R1=irf_resamp(R1,B1);
R2=irf_resamp(R2,B1);
R3=irf_resamp(R3,B1);
R4=irf_resamp(R4,B1);

gradB=c_4_grad('R?','B?','grad');
disp('Calculating null position');
%Calculate the null position
dR1=zeros(length(gradB(:,1)),4);
%First need to calculate how far from each satellite the null is using
%taylor expansion
for i=1:length(gradB(:,1))
    if all(isnan(gradB(i,2:end)))
        continue
    end
    deltaB_null=reshape(gradB(i,2:end),3,3);
    dR1(i,2:4)=B1(i,2:4)/(deltaB_null');  %Row calculations
end
%Add time
Time=B1(:,1);
dR1(:,1)=Time;

%Spacecraft box
minX = min(([R1(:,2) R2(:,2) R3(:,2) R4(:,2)]),[],2);
maxX = max(([R1(:,2) R2(:,2) R3(:,2) R4(:,2)]),[],2);
minY = min(([R1(:,3) R2(:,3) R3(:,3) R4(:,3)]),[],2);
maxY = max(([R1(:,3) R2(:,3) R3(:,3) R4(:,3)]),[],2);
minZ = min(([R1(:,4) R2(:,4) R3(:,4) R4(:,4)]),[],2);
maxZ = max(([R1(:,4) R2(:,4) R3(:,4) R4(:,4)]),[],2);

%Position of the null
Rn1 = irf_add(1,R1,-1,dR1);  %R1-dR1


%This sorting ensures the null is within the four satellites because we
%assumed that the Taylor expansion is only valid for a length of d_i, which
%is the maxmimum separation of the satellites.
disp('Sorting based on the null located within the s/c box thus having a maximum length less than one ion inertial length');
sortNullDx=Rn1(:,2) >= minX & Rn1(:,2) <= maxX;
sortNullDy=Rn1(:,3) >= minY & Rn1(:,3) <= maxY;
sortNullDz=Rn1(:,4) >= minZ & Rn1(:,4) <= maxZ;
sortdr=sortNullDx & sortNullDy & sortNullDz;
%Removing data that didn't fullfill the requirement that the null is
%within ion inertial length distance

if all(~sortdr)   %If none of the sort is fullfilled then just continue
    disp('no null located within the s/c tetrahedron')
    continue;
end

%B
disp('Nulls found within the s/c tetrahedron for time interval')
timestep=median(diff(B1(:,1)));
B1(~sortdr,2:4)=NaN;
tt=B1(sortdr,1);
indStart=tt(1);
indEnd=tt(end);
if indStart==indEnd
    indStart=tt(1)-timestep/2;
    indEnd=tt(end)+timestep/2;
end
disp('Saving time interval')
timeIntervalWithinInterval(iTimeInterval2,:)=[indStart-timestep/2 indEnd+timestep/2];
% Checks if how much of the time interval that didn't fulfill the
% requirement. If it is less than 1/3 of the time interval then the distance to the null is
% assumed to be acceptable
if sum(isnan(B1(:,2))) <= ((1/3)*length(B1(:,1)))
nullPositionFulfilled(iTimeInterval2,1)=true;
end
end

timeIntervalsOfInterest(~nullPositionFulfilled,:)          = NaN;
timeIntervalWithinInterval(timeIntervalWithinInterval==0)  = NaN;
save('timetablenull2003.mat','timeIntervalsOfInterest','timeIntervalWithinInterval');
end