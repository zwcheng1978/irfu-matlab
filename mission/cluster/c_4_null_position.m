function [constraint,Null_Position, c_4_limits,dRmin,Null_Type]=c_4_null_position(B1,B2,B3,B4,R1,R2,R3,R4,varargin)
%C_4_NULL_POSITION - Calculates the null position using 4 spacecraft technique
%
%This function calculates the null position by assuming that the B-field in
%the vicinity of the null can be determined by considering a Taylor
%expansion of the lowest order of B about the null.
%
%   [constraint,Null_Position, c_4_limits,dRmin,Null_Type]=c_4_null_position(B1,B2,B3,B4,R1,R2,R3,R4);
%   [constraint,Null_Position, c_4_limits,dRmin,Null_Type]=c_4_null_position(B1,B2,B3,B4,R1,R2,R3,R4, 'threshold',threshold_value);
%   -threshold=100 means no restriction
%   Null_position = [Time xn yn zn]
%   constraint is a logical vector that shows true when the threshold
%   dRmin = [Time dRmin] Gives the minimum distance to the null point
%   looking from all satellites
%   c_4_limits is a structure containing the maxmimum and minimum positions
%   of all satellites at each time tag.
%   Null_Type is a structure containing each nulltype identified at each
%   time tag.
%   B? is the B-field measured at satellite ?
%   R? is the satellite ? position in GSM coordinate system
%
% Important: Time tags should be included with the vectors and threshold
% should be given in percentage not deciamalform
%
%See Also C_4_NULL_TYPE

%--------written by E.Eriksson--------------------------------------------
n=size(B1,2);
if n<4
    error('Time tag must be included in each input vector. Please do so and try again.')
end
if nargin==0
    help c_4_null_position;
    return;
elseif nargin < 8
    error('Too few input values. See usage: help c_4_null_position')
elseif nargin>10
    error('Too many input values. See usage: help c_4_null_position')
end
if isempty(varargin)==true
    threshold=40;
else
    threshold=varargin{2};
end
%First the magn. field and location of the s/c's need to be
%synchronised and resampled to the same time if there's no NaN values.
if  isempty(isnan(B1(:,1)))
B2=irf_resamp(B2,B1);
B3=irf_resamp(B3,B1);
B4=irf_resamp(B4,B1);
R1=irf_resamp(R1,B1);
R2=irf_resamp(R2,B1);
R3=irf_resamp(R3,B1);
R4=irf_resamp(R4,B1);
end

%Calculates the gradB used in the taylor expansion B0=(B1+B2+B3+B4)/4
gradB=c_4_grad('R?','B?','grad');

%Calculate the null position
dR1=zeros(length(gradB(:,1)),4);
dR2=zeros(length(gradB(:,1)),4);
dR3=zeros(length(gradB(:,1)),4);
dR4=zeros(length(gradB(:,1)),4);
%First need to calculate how far from each satellite the null is using Newton-Rapson method
for i=1:length(gradB(:,1))
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

%The length from each satellite to the null
dRmag1=irf_abs(dR1);
dRmag2=irf_abs(dR2);
dRmag3=irf_abs(dR3);
dRmag4=irf_abs(dR4);

%The minimum distance to the null (which ever satellite that has it)
dRmin(:,2)=min([dRmag1(:,5) dRmag2(:,5) dRmag3(:,5) dRmag4(:,5)], [], 2);
dRmin(:,1)=Time; %adds time

%Null position. They all give the same answer since the null is in one
%position. Calculating each just to check that they do give the same
%position
Rn_1=irf_add(1,R1,-1,dR1);  %R1-dR1 etc.

%Null coordinates
xn=Rn_1(:,2);
yn=Rn_1(:,3);
zn=Rn_1(:,4);
Rn= [xn yn zn];
%Calculate the minimum and maximum values for all s/c's in each direction
%to see distance between s/c's
%Taylor's expansion unreliable distance more than one ion inertial length
%(\approx 1000km)
minX=min(([R1(:,2) R2(:,2) R3(:,2) R4(:,2)]),[],2);
maxX=max(([R1(:,2) R2(:,2) R3(:,2) R4(:,2)]),[],2);
minY=min(([R1(:,3) R2(:,3) R3(:,3) R4(:,3)]),[],2);
maxY=max(([R1(:,3) R2(:,3) R3(:,3) R4(:,3)]),[],2);
minZ=min(([R1(:,4) R2(:,4) R3(:,4) R4(:,4)]),[],2);
maxZ=max(([R1(:,4) R2(:,4) R3(:,4) R4(:,4)]),[],2);

%Check which type the nulls are
[constraint,Nulls]=c_4_null_type(threshold,B1,B2,B3,B4,R1,R2,R3,R4);
%threshold=0.25;
%[eigA,eigB,eigAs,eigBs,eigo,eigx,unknown]=null_type_threshold(B1,B2,B3,B4,R1,R2,R3,R4,threshold)

%For each eigenvalue corresponding to the tolerance level (the two errors less or equal to 40%) break out their corresponding time and dR value
%(the minimum distance from all s/c to the null)

% min and max for all s/c's
minX(constraint,:)=NaN;
maxX(constraint,:)=NaN;
minY(constraint,:)=NaN;
maxY(constraint,:)=NaN;
minZ(constraint,:)=NaN;
maxZ(constraint,:)=NaN;

%Spacecraft max,min position
c_4_limits.minX=minX;
c_4_limits.maxX=maxX;
c_4_limits.minY=minY;
c_4_limits.maxY=maxY;
c_4_limits.minZ=minZ;
c_4_limits.maxZ=maxZ;

%The null positions
Time=B1(:,1);
Time(constraint,:)=NaN;
Rn(constraint,:)=NaN;
Null_Position=[Time Rn];

%dRmin
dRmin(constraint,:)=NaN;

%A
A_null=dRmin;
A_null(Nulls.eigA,:)=NaN;
%B
B_null=dRmin;
B_null(Nulls.eigB,:)=NaN;
%X
x_null=dRmin;
x_null(Nulls.eigx,:)=NaN;
%Bs
Bs_null=dRmin;
Bs_null(Nulls.eigBs,:)=NaN;
%As
As_null=dRmin;
As_null(Nulls.eigAs,:)=NaN;
%O
o_null=dRmin;
o_null(Nulls.eigo,:)=NaN;
%Unknown type
unknown_null=dRmin;
unknown_null(Nulls.unknown,:)=NaN;

Null_Type.A=A_null;
Null_Type.B=B_null;
Null_Type.As=As_null;
Null_Type.Bs=Bs_null;
Null_Type.unknown=unknown_null;
Null_Type.x=x_null;
Null_Type.o=o_null;
end