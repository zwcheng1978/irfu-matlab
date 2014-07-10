function [timeIntervalWithinInterval, timeIntervalsOfInterest]=c_4_null_sorting(timeIntervalsOfInterest,varargin) 
%C_4_NULL_SORTING - Sorts data for finding nulls in the magnetotail
%
%This function sorts data by only including data measured in the
%magnetotail that fulfills the restrictions of maximum spacecraft
%separation of one ion inertial length and that the position of the null is
%within the spacecraft (maximum distance of one ion inertial length).
%
%   [timeIntervalWithinInterval, timeIntervalsOfInterest]=c_4_NULL_SORTING(timeIntervalsOfInterest)
%   [timeIntervalWithinInterval, timeIntervalsOfInterest]=c_4_NULL_SORTING(timeIntervalsOfInterest,'SPIN')
%   OUTPUT
%   B? = the B values at the times that fulfills the restrictions: column 1 - time
%   column 2-4 B-field in x,y,z direction
%   R? = The spacecraft positions at the times that fulfills the restrictions: column 1 - time
%   column 2-4 satellite position in x,y,z direction
%   INPUT
%   'Re-DO' is what is used when the timeIntervalWithinInterval is re-run
%   in the program as the timeIntervalOfInterest. The s/c and B-field
%   with full resolution is used in this phase
%   'SPIN' calls on the Bfield measurement with spin resolution otherwise
%   5VPS resolution is used
%
% Important: Time tags should be included with the vectors and threshold
% should be given in percentage not deciamalform
%
%See Also C_4_NULL_POSITION
%

%Removes time intervals (tint) where the separation between the s/c's weren't
%acceptable in sc_separation_sorting function

%Sorting data so that the maximum distance to a null is smaller or equal to the
%largest distance between all satellites thus keeping it to an ion
%inertial length where Taylor expansion is assumed to work (Bfield behaves
%linearly)

disp('Starting calculation based on the magnetic field and position of the null');
nullPositionFulfilled = false(length(timeIntervalsOfInterest(:,1)),1);
timeIntervalWithinInterval = zeros(length(timeIntervalsOfInterest(:,1)),2);
for iTimeInterval2=1:length(timeIntervalsOfInterest(:,1))
tint=timeIntervalsOfInterest(iTimeInterval2,:); %Goes through each tint every loop
iTimeInterval2
if all(isnan(tint)) %If the tint has been NaN (the separation between the 
    %satellites wasn't fulfilled) continue to next tint
    continue
end
if isempty(varargin)
B1=local.c_read('B_vec_xyz_gse__C1_CP_FGM_5VPS',tint);
B1 = irf_gse2gsm(B1); %Gives the magnetic field of Cluster 1 in GSM
B2=local.c_read('B_vec_xyz_gse__C2_CP_FGM_5VPS',tint);
B2 = irf_gse2gsm(B2); %Gives the magnetic field of Cluster 1 in GSM
B3=local.c_read('B_vec_xyz_gse__C3_CP_FGM_5VPS',tint);
B3 = irf_gse2gsm(B3); %Gives the magnetic field of Cluster 1 in GSM
B4=local.c_read('B_vec_xyz_gse__C4_CP_FGM_5VPS',tint);
B4 = irf_gse2gsm(B4); %Gives the magnetic field of Cluster 1 in GSM
R1 = local.c_read('R1',tint);
R2 = local.c_read('R2',tint);
R3 = local.c_read('R3',tint);
R4 = local.c_read('R4',tint);
R1=irf_gse2gsm(R1);
R2=irf_gse2gsm(R2);
R3=irf_gse2gsm(R3);
R4=irf_gse2gsm(R4);
elseif strcmp(varargin{1}, 'SPIN')
B1=local.c_read('B_vec_xyz_gse__C1_CP_FGM_SPIN',tint);
B1 = irf_gse2gsm(B1); %Gives the magnetic field of Cluster 1 in GSM
B2=local.c_read('B_vec_xyz_gse__C2_CP_FGM_SPIN',tint);
B2 = irf_gse2gsm(B2); %Gives the magnetic field of Cluster 1 in GSM
B3=local.c_read('B_vec_xyz_gse__C3_CP_FGM_SPIN',tint);
B3 = irf_gse2gsm(B3); %Gives the magnetic field of Cluster 1 in GSM
B4=local.c_read('B_vec_xyz_gse__C4_CP_FGM_SPIN',tint);
B4 = irf_gse2gsm(B4); %Gives the magnetic field of Cluster 1 in GSM
R1 = local.c_read('R1',tint);
R2 = local.c_read('R2',tint);
R3 = local.c_read('R3',tint);
R4 = local.c_read('R4',tint);
R1=irf_gse2gsm(R1);
R2=irf_gse2gsm(R2);
R3=irf_gse2gsm(R3);
R4=irf_gse2gsm(R4);
elseif strcmp(varargin{1}, 'RE-DO')
B1=local.c_read('B_vec_xyz_gse__C1_CP_FGM_FULL',tint);
B1 = irf_gse2gsm(B1); %Gives the magnetic field of Cluster 1 in GSM
B2=local.c_read('B_vec_xyz_gse__C2_CP_FGM_FULL',tint);
B2 = irf_gse2gsm(B2); %Gives the magnetic field of Cluster 1 in GSM
B3=local.c_read('B_vec_xyz_gse__C3_CP_FGM_FULL',tint);
B3 = irf_gse2gsm(B3); %Gives the magnetic field of Cluster 1 in GSM
B4=local.c_read('B_vec_xyz_gse__C4_CP_FGM_FULL',tint);
B4 = irf_gse2gsm(B4); %Gives the magnetic field of Cluster 1 in GSM
R1 = local.c_read('sc_pos_xyz_gse__C1_CP_FGM_FULL',tint);
R2 = local.c_read('sc_pos_xyz_gse__C2_CP_FGM_FULL',tint);
R3 = local.c_read('sc_pos_xyz_gse__C3_CP_FGM_FULL',tint);
R4 = local.c_read('sc_pos_xyz_gse__C4_CP_FGM_FULL',tint);
R1=irf_gse2gsm(R1);
R2=irf_gse2gsm(R2);
R3=irf_gse2gsm(R3);
R4=irf_gse2gsm(R4);
end

if isempty(B1)
    tint
    continue
end
if isempty(B2)
    tint
    continue
end
if isempty(B3)
    tint
    continue
end
if isempty(B4)
    tint
    continue
end
%Resample B- and R-fields
B2=irf_resamp(B2,B1);
B3=irf_resamp(B3,B1);
B4=irf_resamp(B4,B1);
R1=irf_resamp(R1,B1);
R2=irf_resamp(R2,B1);
R3=irf_resamp(R3,B1);
R4=irf_resamp(R4,B1);

%Bx=[B1(:,2) B2(:,2) B3(:,2) B4(:,2)];
%By=[B1(:,3) B2(:,3) B3(:,3) B4(:,3)];
%Bz=[B1(:,4) B2(:,4) B3(:,4) B4(:,4)];
%signOfBxfields = sign(Bx);
%allNegBxfieldsHaveSameSign = all(signOfBxfields==-1,2);
%allPosBxfieldsHaveSameSign = all(signOfBxfields==1,2);
%signOfByfields = sign(By);
%allNegByfieldsHaveSameSign = all(signOfByfields==-1,2);
%allPosByfieldsHaveSameSign = all(signOfByfields==1,2);
%allNegBzfieldsHaveSameSign = all(signOfBzfields==-1,2);
%allPosBzfieldsHaveSameSign = all(signOfBzfields==1,2);
%allBfieldsNotAroundZero = allNegBxfieldsHaveSameSign | allPosBxfieldsHaveSameSign ...
%| allNegByfieldsHaveSameSign| allPosByfieldsHaveSameSign | ...
%
%B1(allBfieldsNotAroundZero,2:4)=NaN;
%B2(allBfieldsNotAroundZero,2:4)=NaN;
%B3(allBfieldsNotAroundZero,2:4)=NaN;
%B4(allBfieldsNotAroundZero,2:4)=NaN;

%if sum(isnan(B1(:,2)))==length(B1(:,2)) %If all B values aren't around zero then 
    %no point in trying to find a null
    %continue
%end

%R1(allBfieldsNotAroundZero,2:4)=NaN;
%R2(allBfieldsNotAroundZero,2:4)=NaN;
%R3(allBfieldsNotAroundZero,2:4)=NaN;
%R4(allBfieldsNotAroundZero,2:4)=NaN;

gradB=c_4_grad('R?','B?','grad');
disp('Calculating null position');
%Calculate the null position
dR1=zeros(length(gradB(:,1)),4);
dR2=zeros(length(gradB(:,1)),4);
dR3=zeros(length(gradB(:,1)),4);
dR4=zeros(length(gradB(:,1)),4);
%First need to calculate how far from each satellite the null is using
%taylor expansion
for i=1:length(gradB(:,1))
    if all(isnan(gradB(i,2:end)))
        continue
    end
    deltaB_null=reshape(gradB(i,2:end),3,3);
    dR1(i,2:4)=B1(i,2:4)/(deltaB_null');  %Row calculations
    dR2(i,2:4)=B2(i,2:4)/(deltaB_null'); 
    dR3(i,2:4)=B3(i,2:4)/(deltaB_null'); 
    dR4(i,2:4)=B4(i,2:4)/(deltaB_null'); 
end
%Add time
Time=B1(:,1);
dR1(:,1)=Time;
dR2(:,1)=Time;
dR3(:,1)=Time;
dR4(:,1)=Time;
 %Position of Null
%The length from each satellite to the null
dRmag1 = irf_abs(dR1);
dRmag2 = irf_abs(dR2);
dRmag3 = irf_abs(dR3);
dRmag4 = irf_abs(dR4);

dRmin=zeros(length(dRmag1(:,1)),2);
%The minimum distance to the null (which ever satellite that has it)
dRmin(:,2) = min([dRmag1(:,5) dRmag2(:,5) dRmag3(:,5) dRmag4(:,5)], [], 2);
dRmin(:,1) = Time;

minX = min(([R1(:,2) R2(:,2) R3(:,2) R4(:,2)]),[],2);
maxX = max(([R1(:,2) R2(:,2) R3(:,2) R4(:,2)]),[],2);
minY = min(([R1(:,3) R2(:,3) R3(:,3) R4(:,3)]),[],2);
maxY = max(([R1(:,3) R2(:,3) R3(:,3) R4(:,3)]),[],2);
minZ = min(([R1(:,4) R2(:,4) R3(:,4) R4(:,4)]),[],2);
maxZ = max(([R1(:,4) R2(:,4) R3(:,4) R4(:,4)]),[],2);

Rn1 = irf_add(1,R1,-1,dR1);  %R1-dR1


%This sorting ensures the null is within the four satellites because we
%assumed that the Taylor expansion is only valid for a length of d_i, which
%is the maxmimum separation of the satellites.
d_i=1000; %Ion intertial length
disp('Sorting based on the null located within the s/c tetrahedron thus having a maximum length of one ion inertial length');
sortNull=dRmin(:,2) <= d_i;
sortNullDx=Rn1(:,2) >= minX & Rn1(:,2) <= maxX;
sortNullDy=Rn1(:,3) >= minY & Rn1(:,3) <= maxY;
sortNullDz=Rn1(:,4) >= minZ & Rn1(:,4) <= maxZ;
sortdr= sortNull & sortNullDx & sortNullDy & sortNullDz;
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
    ttindex=iTimeInterval2
nullPositionFulfilled(iTimeInterval2,1)=true;
end
end

timeIntervalsOfInterest(~nullPositionFulfilled,:)          = NaN;
timeIntervalWithinInterval(timeIntervalWithinInterval==0)  = NaN;
save('timetablenull.mat','timeIntervalsOfInterest','timeIntervalWithinInterval');
end