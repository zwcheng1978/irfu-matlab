function [xn,yn,zn, minX, maxX, minY, maxY, minZ, maxZ,dRmin,A_null,B_null,As_null,Bs_null,x_null,o_null,unknown_null]=c_4_null_position(B1,B2,B3,B4,R1,R2,R3,R4)
%c_4_null_position - Calculates the null position using 4 spacecraft technique
%
%This function calculates the null position by assuming that the B-field in
%the vicinity of the null can be determined by considering a Taylor
%expansion of the lowest order of B about the null.
%
%   [xn,yn,zn]=c_4_null_position(B1,B2,B3,B4,R1,R2,R3,R4)
%
% Important: Time tags should be included with the vectors
%
%See Also POINCARE_INDEX_3D, PLOT_NULL2

%--------written by E.Eriksson--------------------------------------------
n=size(B1,2); 
if n<4
error('Time tag must be included in each input vector. Please do so and try again.')
elseif nargin<8
    error('Too few input values. See usage:')
    help null_eig;
    elseif nargin>8
    error('Too many input values. See usage:')
    help null_eig;
else
    %First the magn. field and location of the s/c's need to be
    %synchronised and resampled to the same time.
B2=irf_resamp(B2,B1);
B3=irf_resamp(B3,B1);
B4=irf_resamp(B4,B1);
R1=irf_resamp(R1,B1);
R2=irf_resamp(R2,B1);
R3=irf_resamp(R3,B1);
R4=irf_resamp(R4,B1);

%Calculates the gradB used in the taylor expansion B0=(B1+B2+B3+B4)/4
gradB=c_4_grad('R?','B?','grad');

%Calculate the null position

%First need to calculate how far from each satellite the null is using Newton-Rapson method 
for i=1:length(gradB(:,1))
    deltaB_null=reshape(gradB(i,2:end),3,3);
    dR1(i,2:4)=B1(i,2:4)*inv(deltaB_null');  %Row calculations   
    dR2(i,2:4)=B2(i,2:4)*inv(deltaB_null'); 
    dR3(i,2:4)=B3(i,2:4)*inv(deltaB_null'); 
    dR4(i,2:4)=B4(i,2:4)*inv(deltaB_null'); 
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
Rn_2=irf_add(1,R2,-1,dR2);
Rn_3=irf_add(1,R3,-1,dR3);
Rn_4=irf_add(1,R4,-1,dR4);
Time=B1(:,1);

%Null coordinates
xn=[Time Rn_1(:,2)];
yn=[Time Rn_1(:,3)];
zn=[Time Rn_1(:,4)];

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
[eigA,eigB,eigAs,eigBs,eigo,eigx,unknown]=null_type_ideal(B1,B2,B3,B4,R1,R2,R3,R4);
threshold=0.25;
%[eigA,eigB,eigAs,eigBs,eigo,eigx,unknown]=null_type_threshold(B1,B2,B3,B4,R1,R2,R3,R4,threshold)

%For each eigenvalue break out their corresponding time and dR value
    %(the distance from s/c to the null
    %A
    t_A=dRmin(eigA,1);
    A_dR=dRmin(eigA,2);
    A_null=[t_A A_dR];
    %B
    t_B=dRmin(eigB,1);
    B_dR=dRmin(eigB,2);
    B_null=[t_B B_dR];
    %X
    t_x=dRmin(eigx,1);
    x_dR=dRmin(eigx,2);
    x_null=[t_x x_dR];
    %Bs
    t_Bs=dRmin(eigBs,1);
    Bs_dR=dRmin(eigBs,2);
    Bs_null=[t_Bs Bs_dR];
    %As
    t_As=dRmin(eigAs,1);
    As_dR=dRmin(eigAs,2);
    As_null=[t_As As_dR];
    %O
    t_o=dRmin(eigo,1);
    o_dR=dRmin(eigo,2);
    o_null=[t_o o_dR];
    %Unknown type
    t_unknown=dRmin(unknown,1);
    unknown_dR=dRmin(unknown,2);
    unknown_null=[t_unknown unknown_dR];

%Error in percentage - estimate if linear interpolation is valid to use
[j,divB,B,jxB,divTshear,divPb] = c_4_j('R?','B?');
magn_current=irf_abs(j);
jmag=magn_current(:,[1 5]); %fifth column contains the abs value of j (sqrt(j(:,2).^2+j(:,3).^2+j(:,4).^2)) for each time tag
err_4C=irf_multiply(1,divB,1,jmag,-1); %Essentially does divB/jmag
err_4C(:,2)=abs(err_4C(:,2))*100;

end