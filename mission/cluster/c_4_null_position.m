function [constraint,nullPosition,C4limits,dRmin,dRmax,NullTypeMaxDistance,NullTypeMinDistance]=c_4_null_position(R1,R2,R3,R4,B1,B2,B3,B4,varargin)
%C_4_NULL_POSITION - Calculates the null position using 4 spacecraft technique
%
%This function calculates the null position by assuming that the B-field in
%the vicinity of the null can be determined by considering a Taylor
%expansion of the lowest order of B about the null.
%
%   [constraint,nullPosition,C4limits,dRmin,dRmax,NullTypeMaxDistance,NullTypeMinDistance]...
%   =C_4_NULL_POSITION(R1,R2,R3,R4,B1,B2,B3,B4);
%   [constraint,nullPosition,C4limits,dRmin,dRmax,NullTypeMaxDistance,NullTypeMinDistance]...
%   =C_4_NULL_POSITION(R1,R2,R3,R4,B1,B2,B3,B4, 'threshold',threshold_value);
%   -threshold=100 means no restriction
%   OUTPUT
%   nullPosition = [Time xn yn zn]
%   constraint is a logical vector that shows true when the threshold
%   dRmin = [Time dRmin] Gives the minimum distance to the null point
%   looking from all satellites
%   dRmax = [Time dRmax] Gives the maximum distance to the null point
%   looking from all satellites
%   C4limits is a structure containing the maxmimum and minimum positions
%   of all satellites at each time tag.
%   NullTypeMaxDistance is a structure containing each nulltype identified at each
%   time tag for the maximum distance from the satellites to the nullpoint.
%   NullTypeMinDistance is a structure containing each nulltype identified at each
%   time tag for the minimum distance from the satellites to the nullpoint.
%   INPUT
%   B? = the B-field measured at satellite ?: column 1 - time
%   column 2-4 B-field in x,y,z direction
%   R? = the satellite ? position in GSM coordinate system: column 1 - time
%   column 2-4 satellite position in x,y,z direction
%
% Important: Time tags should be included with the vectors and threshold
% should be given in percentage not deciamalform
%
%See Also C_4_NULL_TYPE

%--------written by E.Eriksson--------------------------------------------
n=size(B1,2);
if n < 4
    error('Time tag must be included in each input vector. Please do so and try again.')
end
if nargin == 0
    help c_4_null_position;
    return;
elseif nargin < 8
    error('Too few input values. See usage: help c_4_null_position')
elseif nargin > 10
    error('Too many input values. See usage: help c_4_null_position')
end
if isempty(varargin) == true
    threshold = 40;
else
    threshold = varargin{2};
end
%First the magn. field and location of the s/c's need to be
%synchronised and resampled to the same time
%if  isempty(isnan(B1(:,1)))
B2 = irf_resamp(B2,B1);
B3 = irf_resamp(B3,B1);
B4 = irf_resamp(B4,B1);
R1 = irf_resamp(R1,B1);
R2 = irf_resamp(R2,B1);
R3 = irf_resamp(R3,B1);
R4 = irf_resamp(R4,B1);
%end

%Calculates the gradB used in the taylor expansion
gradB = c_4_grad('R?','B?','grad');

%Calculate the null position
dR1=zeros(length(gradB(:,1)),4);
dR2=zeros(length(gradB(:,1)),4);
dR3=zeros(length(gradB(:,1)),4);
dR4=zeros(length(gradB(:,1)),4);
%First need to calculate how far from each satellite the null is using
%Taylor Expansion
for i=1:length(gradB(:,1))
    if all(isnan(gradB(i,2:end)))
        continue
    end
    deltaBnull=reshape(gradB(i,2:end),3,3);
    dR1(i,2:4) = B1(i,2:4)/(deltaBnull');  %Row calculations
    dR2(i,2:4) = B2(i,2:4)/(deltaBnull');
    dR3(i,2:4) = B3(i,2:4)/(deltaBnull');
    dR4(i,2:4) = B4(i,2:4)/(deltaBnull');
end
%Add time
Time=B1(:,1);
dR1(:,1)=Time;
dR2(:,1)=Time;
dR3(:,1)=Time;
dR4(:,1)=Time;

%The length from each satellite to the null
dRmag1 = irf_abs(dR1);
dRmag2 = irf_abs(dR2);
dRmag3 = irf_abs(dR3);
dRmag4 = irf_abs(dR4);

%The minimum distance to the null (which ever satellite that has it)
dRmin(:,2) = min([dRmag1(:,5) dRmag2(:,5) dRmag3(:,5) dRmag4(:,5)], [], 2);
dRmin(:,1) = Time; %adds the time

dRmax(:,2) = max([dRmag1(:,5) dRmag2(:,5) dRmag3(:,5) dRmag4(:,5)], [], 2);
dRmax(:,1) = Time; %adds the time

%Null position
Rn1 = irf_add(1,R1,-1,dR1);  %R1-dR1

%Null coordinates
xn = Rn1(:,2);
yn = Rn1(:,3);
zn = Rn1(:,4);
Rn = [xn yn zn];
%Calculate the minimum and maximum values for all s/c's in each direction
%to see distance between s/c's
%Taylor's expansion unreliable distance more than one ion inertial length
%(\approx 1000km)
minX = min(([R1(:,2) R2(:,2) R3(:,2) R4(:,2)]),[],2);
maxX = max(([R1(:,2) R2(:,2) R3(:,2) R4(:,2)]),[],2);
minY = min(([R1(:,3) R2(:,3) R3(:,3) R4(:,3)]),[],2);
maxY = max(([R1(:,3) R2(:,3) R3(:,3) R4(:,3)]),[],2);
minZ = min(([R1(:,4) R2(:,4) R3(:,4) R4(:,4)]),[],2);
maxZ = max(([R1(:,4) R2(:,4) R3(:,4) R4(:,4)]),[],2);

%Check which type the nulls are
[constraint,Nulls] = c_4_null_type(R1,R2,R3,R4,B1,B2,B3,B4,threshold);
%threshold=0.25;
%[eigA,eigB,eigAs,eigBs,eigo,eigx,unknown]=null_type_threshold(B1,B2,B3,B4,R1,R2,R3,R4,threshold)

%For each eigenvalue corresponding to the tolerance level (the two errors less or equal to 40%) break out their corresponding time and dR value
%(the minimum distance from all s/c to the null)

% min and max for all s/c's
minX(~constraint,:) = NaN;
maxX(~constraint,:) = NaN;
minY(~constraint,:) = NaN;
maxY(~constraint,:) = NaN;
minZ(~constraint,:) = NaN;
maxZ(~constraint,:) = NaN;

%Spacecraft max,min position
C4limits.minX = minX;
C4limits.maxX = maxX;
C4limits.minY = minY;
C4limits.maxY = maxY;
C4limits.minZ = minZ;
C4limits.maxZ = maxZ;

%The null positions
Time=B1(:,1);
Rn(~constraint,:) = NaN;
nullPosition      = [Time Rn];

%dRmin
dRmin(~constraint,2)= NaN;

%A
distanceANull                         = dRmin;
distanceANull(~Nulls.eigA,2)          = NaN;
%B
distanceBNull                         = dRmin;
distanceBNull(~Nulls.eigB,2)          = NaN;
%X
distanceXNull                         = dRmin;
distanceXNull(~Nulls.eigx,2)          = NaN;
%Bs
distanceBsNull                        = dRmin;
distanceBsNull(~Nulls.eigBs,2)        = NaN;
%As
distanceAsNull                        = dRmin;
distanceAsNull(~Nulls.eigAs,2)        = NaN;
%O
distanceONull                         = dRmin;
distanceONull(~Nulls.eigo,2)          = NaN;
%Unknown type
distanceUnknownNull                   = dRmin;
distanceUnknownNull(~Nulls.unknown,2) = NaN;

NullTypeMinDistance.A       = distanceANull;
NullTypeMinDistance.B       = distanceBNull;
NullTypeMinDistance.As      = distanceAsNull;
NullTypeMinDistance.Bs      = distanceBsNull;
NullTypeMinDistance.unknown = distanceUnknownNull;
NullTypeMinDistance.x       = distanceXNull;
NullTypeMinDistance.o       = distanceONull;

%dRmax
dRmax(~constraint,2)= NaN;

%A
distanceANull                         = dRmax;
distanceANull(~Nulls.eigA,2)          = NaN;
%B
distanceBNull                         = dRmax;
distanceBNull(~Nulls.eigB,2)          = NaN;
%X
distanceXNull                         = dRmax;
distanceXNull(~Nulls.eigx,2)          = NaN;
%Bs
distanceBsNull                        = dRmax;
distanceBsNull(~Nulls.eigBs,2)        = NaN;
%As
distanceAsNull                        = dRmax;
distanceAsNull(~Nulls.eigAs,2)        = NaN;
%O
distanceONull                         = dRmax;
distanceONull(~Nulls.eigo,2)          = NaN;
%Unknown type
distanceUnknownNull                   = dRmax;
distanceUnknownNull(~Nulls.unknown,2) = NaN;

NullTypeMaxDistance.A       = distanceANull;
NullTypeMaxDistance.B       = distanceBNull;
NullTypeMaxDistance.As      = distanceAsNull;
NullTypeMaxDistance.Bs      = distanceBsNull;
NullTypeMaxDistance.unknown = distanceUnknownNull;
NullTypeMaxDistance.x       = distanceXNull;
NullTypeMaxDistance.o       = distanceONull;
end