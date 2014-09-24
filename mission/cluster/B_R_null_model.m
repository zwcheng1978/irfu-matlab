function [Modeldata]=B_R_null_model(s,p,j)
%B_R_null_model - Creates a simple model of B-field with no currents and position (R)
%measurements for 4 s/c's. 

% This function creates a simple linear Bfield that can have either complex
% or real eigenvalues. The s is the scalar parameter which means the
% eigenvalue is either A-type (s<0) or B-type (s>0). p is a freely chosen
% parameter which gives different eigenvalues. The higher p chosen the
% larger threshold current you have and the higher current in the spine (j)
% you need for As or Bs nulls. Default value is 0.2. j=0 gives the old
% B-model were you only have A and B type nulls. First loop changes s
% values (if you want), second loops over the third loop a 1000 times to see 
%the random functions distribution. The third loop gradually increases the amplitude of
% the random vectorised disturbance in B3.
%This function also saves data from each loop run so a lot of code is just
%the save function.

% Important: Time tags should be included with the vectors
%
%

%Calculates the s/c's position
% % S/C 1-4 space coordinates using an estimated spacing between satellites
% from Cluster in Oct 1 2001 
%Keeping the s/c's stuck and moving the null instead for simplicity sake
if nargin == 0
    help B_R_null_model;
    return;
elseif nargin < 3
    error('Too few input values. See usage: help B_R_null_model')
elseif nargin > 3
    error('Too many input values. See usage: help B_R_null_model')
end

%Creates time axis
t=(0.04:0.04:20)';
t1=1001929700.10000-0.0400;
tt=zeros(length(t),1);
x1=zeros(length(t),1);
x2=zeros(length(t),1);
x3=zeros(length(t),1);
x4=zeros(length(t),1);
y1=zeros(length(t),1);
y2=zeros(length(t),1);
y3=zeros(length(t),1);
y4=zeros(length(t),1);
z1=zeros(length(t),1);
z2=zeros(length(t),1);
z3=zeros(length(t),1);
z4=zeros(length(t),1);
%For each time gives the satellite position which is fixed in time
for i=1:length(t)
    t2=t1+0.04;
    x1(i,1)=-3581.8007812500;
    x2(i,1)=-3342.0477343544;
    x3(i,1)=-3386.8466796677;
    x4(i,1)=-3450.1474609177;
    y1(i,1)=1341.4572387373;
    y2(i,1)=1198.3072639669;
    y3(i,1)=1415.9244835693;
    y4(i,1)=1142.4099062605;
    z1(i,1)=1759.6076039468;
    z2(i,1)=1830.2748247349;
    z3(i,1)=1637.9525007500;
    z4(i,1)=1626.6903388264;
    tt(i,1)=t1+0.0400;
    t1=t2;
end
% S/C 1-4 space coordinates using an estimated spacing between satellites
% from Cluster in Oct 1 2001 
R1=[tt x1 y1 z1];
R2=[tt x2 y2 z2];
R3=[tt x3 y3 z3];
R4=[tt x4 y4 z4];

%Calculates the Bfield measured by S/C's for each cordinate the magnetic
%null is in position (xn,yn,zn)

%% Assumption of a,b,c,alpha and beta (--+ A, ++- B) A +1 B -1
%Moving null
x0=-3200;
y0=1100;
z0=1500;
xn=x0-33.*t;
yn=y0+25.*t;
zn=z0+25.*t;
Rn=[tt xn yn zn];
%ii loops over changes in the scalar parameter (if you want to increase or
%decrease eigenvalues
for ii=1:1%21
    %Threshold current used in B-model
    jthresh=sqrt((p-1).^2);
    %Eigenvalues
    landa1=(p+1+sqrt(jthresh.^2-j.^2))/(s*2);
    landa2=(p+1-sqrt(jthresh.^2-j.^2))/(s*2);
    landa3=-(p+1)/s;
    %iii for loop runs to checks the distribution of the random function
    for iii=1:1000
        %Starting value for the amplitude of the disturbance
    amplitude=0;
    uncertaincyamplitude1=0;
    uncertaincyamplitude2=0;
    uncertaincyamplitude4=0;
    %i for loop loops over each amplitude increase so what i and iii does
    %toghether is calculate for each amplitude the distribution of the
    %random function since it will be changed for each iii run.
    for i=1:500
    %The spherical coordinate angles
    phi=rand(size(tt)).*180;
    theta=rand(size(tt)).*360;
    
    phi1=rand(size(tt)).*180;
    theta1=rand(size(tt)).*360;
    
    phi2=rand(size(tt)).*180;
    theta2=rand(size(tt)).*360;
    
    phi4=rand(size(tt)).*180;
    theta4=rand(size(tt)).*360;
    
    %Disturbance in B-field for each satellite in this run only B3
    %experiences a disturbance
    Bz = amplitude.*cosd(phi);
    By = amplitude.*sind(theta).*sind(phi);
    Bx = amplitude.*cosd(theta).*sind(phi);
    
    Buncertaintyz1 = uncertaincyamplitude1.*cosd(phi1);
    Buncertaintyy1 = uncertaincyamplitude1.*sind(theta1).*sind(phi1);
    Buncertaintyx1 = uncertaincyamplitude1.*cosd(theta1).*sind(phi1);
    
    Buncertaintyz2 = uncertaincyamplitude2.*cosd(phi2);
    Buncertaintyy2 = uncertaincyamplitude2.*sind(theta2).*sind(phi2);
    Buncertaintyx2 = uncertaincyamplitude2.*cosd(theta2).*sind(phi2);
    
    Buncertaintyz4 = uncertaincyamplitude4.*cosd(phi4);
    Buncertaintyy4 = uncertaincyamplitude4.*sind(theta4).*sind(phi4);
    Buncertaintyx4 = uncertaincyamplitude4.*cosd(theta4).*sind(phi4);
    
    %Magnetic field measured at the third satellite
    B3x=(1/s)*(x3-xn)-(1/2)*(1/s)*j.*(y3-yn);
    B3y=(p/s).*(y3-yn)+(1/(s*2))*j.*(x3-xn);
    B3z=-((p+1)/s).*(z3-zn);
    B3x=B3x+Bx;
    B3y=B3y+By;
    B3z=B3z+Bz;
    B3=[tt B3x B3y B3z];
    
    %Magnetic field measured at the first satellite
    B1x=(1/s)*(x1-xn)-(1/2)*(1/s)*j.*(y1-yn);
    B1y=(p/s).*(y1-yn)+(1/(s*2))*j.*(x1-xn);
    B1z=-((p+1)/s).*(z1-zn);
    B1x=B1x+Buncertaintyx1;
    B1y=B1y+Buncertaintyy1;
    B1z=B1z+Buncertaintyz1;
    B1=[tt B1x B1y B1z];
    
    %Magnetic field measured at the second satellite
    B2x=(1/s)*(x2-xn)-(1/2)*(1/s)*j.*(y2-yn);
    B2y=(p/s).*(y2-yn)+(1/(s*2))*j.*(x2-xn);
    B2z=-((p+1)/s).*(z2-zn);
    B2x=B2x+Buncertaintyx2;
    B2y=B2y+Buncertaintyy2;
    B2z=B2z+Buncertaintyz2;
    B2=[tt B2x B2y B2z];
    
    %Magnetic field measured at the fourth satellite
    B4x=(1/s)*(x4-xn)-(1/2)*(1/s)*j.*(y4-yn);
    B4y=(p/s).*(y4-yn)+(1/(s*2))*j.*(x4-xn);
    B4z=-((p+1)/s).*(z4-zn);
    B4x=B4x+Buncertaintyx4;
    B4y=B4y+Buncertaintyy4;
    B4z=B4z+Buncertaintyz4;
    B4=[tt B4x B4y B4z];
    
    %Calculates for each time step the size of the B-field difference for
    %all satellites
    mBx=max([B1(:,2) B2(:,2) B3(:,2) B4(:,2)],[],2)-min([B1(:,2) B2(:,2) B3(:,2) B4(:,2)],[],2);
    mBy=max([B1(:,3) B2(:,3) B3(:,3) B4(:,3)],[],2)-min([B1(:,3) B2(:,3) B3(:,3) B4(:,3)],[],2);
    mBz=max([B1(:,4) B2(:,4) B3(:,4) B4(:,4)],[],2)-min([B1(:,4) B2(:,4) B3(:,4) B4(:,4)],[],2);
    maxBx=mBx;
    maxBy=mBy;
    maxBz=mBz;
    Bmag1 = irf_abs(B1);
    Bmag2 = irf_abs(B2);
    Bmag3 = irf_abs(B3);
    Bmag4 = irf_abs(B4);
    %Calculates for each time step the size between the magnitude of all
    %satellites B-fields.
    maxBtot = max([Bmag1(:,5) Bmag2(:,5) Bmag3(:,5) Bmag4(:,5)], [], 2)-min([Bmag1(:,5) Bmag2(:,5) Bmag3(:,5) Bmag4(:,5)], [], 2);
    
    %For each amplitude the ampltide and the ratio between the amplitude of
    %the B3 disturbance divided by each bfield size in each direction (in
    %percentage)
    KvotBetweenamplitudeandBxieldDifference=[amplitude ((amplitude./maxBx))'.*100];
    KvotBetweenamplitudeandByieldDifference=[amplitude  ((amplitude./maxBy))'.*100];
    KvotBetweenamplitudeandBzieldDifference=[amplitude  ((amplitude./maxBz))'.*100];
    KvotBetweenamplitudeandBtotDifference=[amplitude  ((amplitude./maxBtot))'.*100];
    
    %Calculates the nullPosition and the NullType
   [nullPosition,~,~,NullType,~,~,~,~,~]=c_4_null_position(R1,R2,R3,R4,B1,B2,B3,B4); 
   
   %Here starts the saving function
    nullPositionlogical = false(length(nullPosition(:,1)),1);
    
    for TimeInterval2=1:length(nullPosition(:,1))
        if all(isnan(nullPosition(TimeInterval2,2:4))) %If the nullPosition has been NaN continue to the next nullPosition
            continue
        end
        nullPositionlogical(TimeInterval2,1)=true;
    end
    disp('calculates number of type etc.');
   % nullPosition(~nullPositionlogical,:)=NaN;
   %saves only the Rn where we've located a null.
    Rn(~nullPositionlogical,:)=NaN;
    
    %sigmaerror=amplitude.*ones(size(nullPosition(:,2)),1);
    %xerror=[sigmaerror tt Rn(:,2)-nullPosition(:,2)];
    %yerror=[sigmaerror tt Rn(:,3)-nullPosition(:,3)];
    %zerror=[sigmaerror tt Rn(:,4)-nullPosition(:,4)];
    
    %amplitudesave=amplitude.*ones(size(B1(:,2)),1);
    %time=B1(:,1);
    %Bfield1=[time amplitudesave B1(:,2:4)];
    %Bfield2=[time amplitudesave B2(:,2:4)];
    %Bfield3=[time amplitudesave B3(:,2:4)];
    %Bfield4=[time amplitudesave B4(:,2:4)];
    
    %Saves the number of each type found for each amplitude 500 is the size
    %of the Bfields so I remove all NaN values where a certain type hasn't
    %been found.
    A=[amplitude 500-sum(isnan(NullType.position.A(:,2)))];
    B=[amplitude 500-sum(isnan(NullType.position.B(:,2)))];
    As=[amplitude 500-sum(isnan(NullType.position.As(:,2)))];
    Bs=[amplitude 500-sum(isnan(NullType.position.Bs(:,2)))];
    unknown=[amplitude 500-sum(isnan(NullType.position.unknown(:,2)))];
    x=[amplitude 500-sum(isnan(NullType.position.x(:,2)))];
    o=[amplitude 500-sum(isnan(NullType.position.o(:,2)))];
    numberofnullserror=[amplitude A(1,2)+B(1,2)+Bs(1,2)+As(1,2)];
    
    %Increases the amplitude with each i step so that in the end you will
    %have an amplitude which is about 75% of mBx. This is so the amplitude
    %is big enought for each increase of decrease when changing s in iii.
    amplitude=amplitude+((mBx(1,1)*0.75)/500);
    %Not causing disturbances at the other satellites but it's possible to
    %implement.
    %uncertaincyamplitude1=uncertaincyamplitude1+0.005;
    %uncertaincyamplitude2=uncertaincyamplitude1+0.005;
    %uncertaincyamplitude4=uncertaincyamplitude1+0.005;
    
    %Just saving function to save data for each run
    if i==1
    KvotAmplitudeandBfieldDiff.Bx=KvotBetweenamplitudeandBxieldDifference;
    KvotAmplitudeandBfieldDiff.By=KvotBetweenamplitudeandByieldDifference;
    KvotAmplitudeandBfieldDiff.Bz=KvotBetweenamplitudeandBzieldDifference;
    KvotAmplitudeandBfieldDiff.Btot=KvotBetweenamplitudeandBtotDifference;
    
    %Bfield.B1=Bfield1;
    %Bfield.B2=Bfield2;
    %Bfield.B3=Bfield3;
    %Bfield.B4=Bfield4;
    
    NumberofType.A=A;
    NumberofType.As=As;
    NumberofType.B=B;
    NumberofType.Bs=Bs;
    NumberofType.unknown=unknown;
    NumberofType.x=x;
    NumberofType.o=o;
    
    NumberofNullsError=numberofnullserror;
    
    %PositionError.x=xerror;
    %PositionError.y=yerror;
    %PositionError.z=zerror;
    
    %Saves data for each next i step in the next row for the different
    %amplitudes
    else
    KvotAmplitudeandBfieldDiff.Bx(length(KvotAmplitudeandBfieldDiff.Bx(:,1))+1:length(KvotBetweenamplitudeandBxieldDifference(:,1))+length(KvotAmplitudeandBfieldDiff.Bx(:,1)),:)=KvotBetweenamplitudeandBxieldDifference;    
    KvotAmplitudeandBfieldDiff.By(length(KvotAmplitudeandBfieldDiff.By(:,1))+1:length(KvotBetweenamplitudeandByieldDifference(:,1))+length(KvotAmplitudeandBfieldDiff.By(:,1)),:)=KvotBetweenamplitudeandByieldDifference;
    KvotAmplitudeandBfieldDiff.Bz(length(KvotAmplitudeandBfieldDiff.Bz(:,1))+1:length(KvotBetweenamplitudeandBzieldDifference(:,1))+length(KvotAmplitudeandBfieldDiff.Bz(:,1)),:)=KvotBetweenamplitudeandBzieldDifference;
    KvotAmplitudeandBfieldDiff.Btot(length(KvotAmplitudeandBfieldDiff.Btot(:,1))+1:length(KvotBetweenamplitudeandBtotDifference(:,1))+length(KvotAmplitudeandBfieldDiff.Btot(:,1)),:)=KvotBetweenamplitudeandBtotDifference;
    
    %Bfield.B1(length(Bfield.B1(:,1))+1:length(Bfield1(:,1))+length(Bfield.B1(:,1)),:)=Bfield1;
    %Bfield.B2(length(Bfield.B2(:,1))+1:length(Bfield2(:,1))+length(Bfield.B2(:,1)),:)=Bfield2;
    %Bfield.B3(length(Bfield.B3(:,1))+1:length(Bfield3(:,1))+length(Bfield.B3(:,1)),:)=Bfield3;
    %Bfield.B4(length(Bfield.B4(:,1))+1:length(Bfield4(:,1))+length(Bfield.B4(:,1)),:)=Bfield4;
    
    NumberofType.A(length(NumberofType.A(:,1))+1:length(A(:,1))+length(NumberofType.A(:,1)),:)=A;
    NumberofType.As(length(NumberofType.As(:,1))+1:length(As(:,1))+length(NumberofType.As(:,1)),:)=As;
    NumberofType.B(length(NumberofType.B(:,1))+1:length(B(:,1))+length(NumberofType.B(:,1)),:)=B;
    NumberofType.Bs(length(NumberofType.Bs(:,1))+1:length(Bs(:,1))+length(NumberofType.Bs(:,1)),:)=Bs;
    NumberofType.x(length(NumberofType.x(:,1))+1:length(x(:,1))+length(NumberofType.x(:,1)),:)=x;
    NumberofType.o(length(NumberofType.o(:,1))+1:length(o(:,1))+length(NumberofType.o(:,1)),:)=o;
    NumberofType.unknown(length(NumberofType.unknown(:,1))+1:length(unknown(:,1))+length(NumberofType.unknown(:,1)),:)=unknown;
    
    NumberofNullsError(length(NumberofNullsError(:,1))+1:length(numberofnullserror(:,1))+length(NumberofNullsError(:,1)),:)=numberofnullserror;
    
    %PositionError.x(length(PositionError.x(:,1))+1:length(xerror(:,1))+length(PositionError.x(:,1)),:)=xerror;
   % PositionError.y(length(PositionError.y(:,1))+1:length(yerror(:,1))+length(PositionError.y(:,1)),:)=yerror;
    %PositionError.z(length(PositionError.z(:,1))+1:length(zerror(:,1))+length(PositionError.z(:,1)),:)=zerror;
    end
    end
     
    %Saving data for each iii run (random function distribution)
if iii==1
    Modeldata.KvotAmplitudeandBfieldDiff.Bx=KvotAmplitudeandBfieldDiff.Bx;
    Modeldata.KvotAmplitudeandBfieldDiff.By=KvotAmplitudeandBfieldDiff.By;
    Modeldata.KvotAmplitudeandBfieldDiff.Bz=KvotAmplitudeandBfieldDiff.Bz;
    Modeldata.KvotAmplitudeandBfieldDiff.Btot=KvotAmplitudeandBfieldDiff.Btot;
    
   % Modeldata.Bfield.B1.Bx=[Bfield.B1(:,1) Bfield.B1(:,2) Bfield.B1(:,3)];
    %Modeldata.Bfield.B1.By=[Bfield.B1(:,1) Bfield.B1(:,2) Bfield.B1(:,4)];
    %Modeldata.Bfield.B1.Bz=[Bfield.B1(:,1) Bfield.B1(:,2) Bfield.B1(:,5)];
    %Modeldata.Bfield.B2.Bx=[Bfield.B2(:,1) Bfield.B2(:,2) Bfield.B2(:,3)];
    %Modeldata.Bfield.B2.By=[Bfield.B2(:,1) Bfield.B2(:,2) Bfield.B2(:,4)];
    %Modeldata.Bfield.B2.Bz=[Bfield.B2(:,1) Bfield.B2(:,2) Bfield.B2(:,5)];
    %Modeldata.Bfield.B3.Bx=[Bfield.B3(:,1) Bfield.B3(:,2) Bfield.B3(:,3)];
    %Modeldata.Bfield.B3.By=[Bfield.B3(:,1) Bfield.B3(:,2) Bfield.B3(:,4)];
    %Modeldata.Bfield.B3.Bz=[Bfield.B3(:,1) Bfield.B3(:,2) Bfield.B3(:,5)];
    %Modeldata.Bfield.B4.Bx=[Bfield.B4(:,1) Bfield.B4(:,2) Bfield.B4(:,3)];
    %Modeldata.Bfield.B4.By=[Bfield.B4(:,1) Bfield.B4(:,2) Bfield.B4(:,4)];
    %Modeldata.Bfield.B4.Bz=[Bfield.B4(:,1) Bfield.B4(:,2) Bfield.B4(:,5)];
    
    Modeldata.NumberofType.A=NumberofType.A;
    Modeldata.NumberofType.As=NumberofType.As;
    Modeldata.NumberofType.B=NumberofType.B;
    Modeldata.NumberofType.Bs=NumberofType.Bs;
    Modeldata.NumberofType.unknown=NumberofType.unknown;
    Modeldata.NumberofType.x=NumberofType.x;
    Modeldata.NumberofType.o=NumberofType.o;
    
    Modeldata.Eigenvalues=[landa1 landa2 landa3];
    
    Modeldata.NumberofNullsError=NumberofNullsError;
    
    Modeldata.scale.s=s;
    Modeldata.scale.p=p;
    Modeldata.scale.j=j;
    
    
    %Modeldata.PositionError.x=PositionError.x;
    %Modeldata.PositionError.y=PositionError.y;
    %Modeldata.PositionError.z=PositionError.z;
     
    %Saves data for each next iii step in the next column for random function
    %distribution
else
    disp('Saves data into Modeldata file');
    Modeldata.KvotAmplitudeandBfieldDiff.Bx(:,length(Modeldata.KvotAmplitudeandBfieldDiff.Bx(1,:))+1:1+length(Modeldata.KvotAmplitudeandBfieldDiff.Bx(1,:)))=KvotAmplitudeandBfieldDiff.Bx(:,2);
    Modeldata.KvotAmplitudeandBfieldDiff.By(:,length(Modeldata.KvotAmplitudeandBfieldDiff.By(1,:))+1:1+length(Modeldata.KvotAmplitudeandBfieldDiff.By(1,:)))=KvotAmplitudeandBfieldDiff.By(:,2);
    Modeldata.KvotAmplitudeandBfieldDiff.Bz(:,length(Modeldata.KvotAmplitudeandBfieldDiff.Bz(1,:))+1:1+length(Modeldata.KvotAmplitudeandBfieldDiff.Bz(1,:)))=KvotAmplitudeandBfieldDiff.Bz(:,2);
    Modeldata.KvotAmplitudeandBfieldDiff.Btot(:,length(Modeldata.KvotAmplitudeandBfieldDiff.Btot(1,:))+1:1+length(Modeldata.KvotAmplitudeandBfieldDiff.Btot(1,:)))=KvotAmplitudeandBfieldDiff.Btot(:,2);
    
   % Modeldata.Bfield.B1.Bx(:,length(Modeldata.Bfield.B1.Bx(1,:))+1:1+length(Modeldata.Bfield.B1.Bx(1,:)))=Bfield.B1(:,3);
    %Modeldata.Bfield.B1.By(:,length(Modeldata.Bfield.B1.By(1,:))+1:1+length(Modeldata.Bfield.B1.By(1,:)))=Bfield.B1(:,4);
    %Modeldata.Bfield.B1.Bz(:,length(Modeldata.Bfield.B1.Bz(1,:))+1:1+length(Modeldata.Bfield.B1.Bz(1,:)))=Bfield.B1(:,5);
    %Modeldata.Bfield.B2.Bx(:,length(Modeldata.Bfield.B2.Bx(1,:))+1:1+length(Modeldata.Bfield.B2.Bx(1,:)))=Bfield.B2(:,3);
    %Modeldata.Bfield.B2.By(:,length(Modeldata.Bfield.B2.By(1,:))+1:1+length(Modeldata.Bfield.B2.By(1,:)))=Bfield.B2(:,4);
   % Modeldata.Bfield.B2.Bz(:,length(Modeldata.Bfield.B2.Bz(1,:))+1:1+length(Modeldata.Bfield.B2.Bz(1,:)))=Bfield.B2(:,5);
    %Modeldata.Bfield.B3.Bx(:,length(Modeldata.Bfield.B3.Bx(1,:))+1:1+length(Modeldata.Bfield.B3.Bx(1,:)))=Bfield.B3(:,3);
    %Modeldata.Bfield.B3.By(:,length(Modeldata.Bfield.B3.By(1,:))+1:1+length(Modeldata.Bfield.B3.By(1,:)))=Bfield.B3(:,4);
    %Modeldata.Bfield.B3.Bz(:,length(Modeldata.Bfield.B3.Bz(1,:))+1:1+length(Modeldata.Bfield.B3.Bz(1,:)))=Bfield.B3(:,5);
    %Modeldata.Bfield.B4.Bx(:,length(Modeldata.Bfield.B4.Bx(1,:))+1:1+length(Modeldata.Bfield.B4.Bx(1,:)))=Bfield.B4(:,3);
   % Modeldata.Bfield.B4.By(:,length(Modeldata.Bfield.B4.By(1,:))+1:1+length(Modeldata.Bfield.B4.By(1,:)))=Bfield.B4(:,4);
   % Modeldata.Bfield.B4.Bz(:,length(Modeldata.Bfield.B4.Bz(1,:))+1:1+length(Modeldata.Bfield.B4.Bz(1,:)))=Bfield.B4(:,5);
    
    Modeldata.NumberofType.A(:,length(Modeldata.NumberofType.A(1,:))+1:1+length(Modeldata.NumberofType.A(1,:)))=NumberofType.A(:,2);
    Modeldata.NumberofType.As(:,length(Modeldata.NumberofType.As(1,:))+1:1+length(Modeldata.NumberofType.As(1,:)))=NumberofType.As(:,2);
    Modeldata.NumberofType.B(:,length(Modeldata.NumberofType.B(1,:))+1:1+length(Modeldata.NumberofType.B(1,:)))=NumberofType.B(:,2);
    Modeldata.NumberofType.Bs(:,length(Modeldata.NumberofType.Bs(1,:))+1:1+length(Modeldata.NumberofType.Bs(1,:)))=NumberofType.Bs(:,2);
    Modeldata.NumberofType.x(:,length(Modeldata.NumberofType.x(1,:))+1:1+length(Modeldata.NumberofType.x(1,:)))=NumberofType.x(:,2);
    Modeldata.NumberofType.o(:,length(Modeldata.NumberofType.o(1,:))+1:1+length(Modeldata.NumberofType.o(1,:)))=NumberofType.o(:,2);
    Modeldata.NumberofType.unknown(:,length(Modeldata.NumberofType.unknown(1,:))+1:1+length(Modeldata.NumberofType.unknown(1,:)))=NumberofType.unknown(:,2);
  
    Modeldata.NumberofNullsError(:,length(Modeldata.NumberofNullsError(1,:))+1:1+length(Modeldata.NumberofNullsError(1,:)))=NumberofNullsError(:,2);
    
    %Modeldata.PositionError.x(:,length(Modeldata.PositionError.x(1,:))+1:1+length(Modeldata.PositionError.x(1,:)))=PositionError.x(:,3);
    %Modeldata.PositionError.y(:,length(Modeldata.PositionError.y(1,:))+1:1+length(Modeldata.PositionError.y(1,:)))=PositionError.y(:,3);
    %Modeldata.PositionError.z(:,length(Modeldata.PositionError.z(1,:))+1:1+length(Modeldata.PositionError.z(1,:)))=PositionError.z(:,3);
end
% Saves data in each run incase something happens
name=['ModelldataA',num2str(ii),'.mat'];
save(name,'Modeldata','-v7.3');
%folder=num2str(ii);
%cd(folder)
%save('ModelldataA.mat','Modeldata','-v7.3');
%cd ..
%iii
    end
    %Saves data for each run with different s values.
if ii==1
    disp('Saves data into Modeldata file');
    ssave=s.*ones(length(Modeldata.KvotAmplitudeandBfieldDiff.Bx(:,1)),1);
    Model.KvotAmplitudeandBfieldDiff.Bx=[ssave Modeldata.KvotAmplitudeandBfieldDiff.Bx];
   Model.KvotAmplitudeandBfieldDiff.By=[ssave Modeldata.KvotAmplitudeandBfieldDiff.By];
   Model.KvotAmplitudeandBfieldDiff.Bz=[ssave Modeldata.KvotAmplitudeandBfieldDiff.Bz];
   Model.KvotAmplitudeandBfieldDiff.Btot=[ssave Modeldata.KvotAmplitudeandBfieldDiff.Btot];
    
   seig=s.*ones(length(Modeldata.Eigenvalues(:,1)),1); 
   Model.Eigenvalues=[seig Modeldata.Eigenvalues];
    
  % sB=s.*ones(length(Modeldata.Bfield.B1.Bx(:,1)),1); 
   % Model.Bfield.B1.Bx=[sB Modeldata.Bfield.B1.Bx];
   % Model.Bfield.B1.By=[sB Modeldata.Bfield.B1.By];
   % Model.Bfield.B1.Bz=[sB Modeldata.Bfield.B1.Bz];
   % Model.Bfield.B2.Bx=[sB Modeldata.Bfield.B2.Bx];
   % Model.Bfield.B2.By=[sB Modeldata.Bfield.B2.By];
   % Model.Bfield.B2.Bz=[sB Modeldata.Bfield.B2.Bz];
   % Model.Bfield.B3.Bx=[sB Modeldata.Bfield.B3.Bx];
   % Model.Bfield.B3.By=[sB Modeldata.Bfield.B3.By];
   % Model.Bfield.B3.Bz=[sB Modeldata.Bfield.B3.Bz];
   % Model.Bfield.B4.Bx=[sB Modeldata.Bfield.B4.Bx];
   % Model.Bfield.B4.By=[sB Modeldata.Bfield.B4.By];
   % Model.Bfield.B4.Bz=[sB Modeldata.Bfield.B4.Bz];
    
    stype=s.*ones(length(Modeldata.NumberofType.A(:,1)),1); 
    Model.NumberofType.A=[stype Modeldata.NumberofType.A];
    Model.NumberofType.As=[stype Modeldata.NumberofType.As];
    Model.NumberofType.B=[stype Modeldata.NumberofType.B];
    Model.NumberofType.Bs=[stype Modeldata.NumberofType.Bs];
    Model.NumberofType.x=[stype Modeldata.NumberofType.x];
    Model.NumberofType.o=[stype Modeldata.NumberofType.o];
    Model.NumberofType.unknown=[stype Modeldata.NumberofType.unknown];
    
    serror=s.*ones(length(Modeldata.NumberofNullsError(:,1)),1); 
    Model.NumberofNullsError=[serror Modeldata.NumberofNullsError];
    
    %sposerrorx=s.*ones(length(Modeldata.PositionError.x(:,1)),1);
    %sposerrory=s.*ones(length(Modeldata.PositionError.y(:,1)),1);
    %sposerrorz=s.*ones(length(Modeldata.PositionError.z(:,1)),1);
    
    %Model.PositionError.x=[sposerrorx Modeldata.PositionError.x];
   % Model.PositionError.y=[sposerrory Modeldata.PositionError.y];
    %Model.PositionError.z=[sposerrorz Modeldata.PositionError.z];
else
     ssave=s.*ones(length(Modeldata.KvotAmplitudeandBfieldDiff.Bx(:,1)),1);
    Model.KvotAmplitudeandBfieldDiff.Bx(1+length(Model.KvotAmplitudeandBfieldDiff.Bx(:,1)):length(Model.KvotAmplitudeandBfieldDiff.Bx(:,1))+length(Modeldata.KvotAmplitudeandBfieldDiff.Bx(:,1)),:)=[ssave Modeldata.KvotAmplitudeandBfieldDiff.Bx];
   Model.KvotAmplitudeandBfieldDiff.By(1+length(Model.KvotAmplitudeandBfieldDiff.By(:,1)):length(Model.KvotAmplitudeandBfieldDiff.By(:,1))+length(Modeldata.KvotAmplitudeandBfieldDiff.By(:,1)),:)=[ssave Modeldata.KvotAmplitudeandBfieldDiff.By];
   Model.KvotAmplitudeandBfieldDiff.Bz(1+length(Model.KvotAmplitudeandBfieldDiff.Bz(:,1)):length(Model.KvotAmplitudeandBfieldDiff.Bz(:,1))+length(Modeldata.KvotAmplitudeandBfieldDiff.Bz(:,1)),:)=[ssave Modeldata.KvotAmplitudeandBfieldDiff.Bz];
  Model.KvotAmplitudeandBfieldDiff.Btot(1+length(Model.KvotAmplitudeandBfieldDiff.Btot(:,1)):length(Model.KvotAmplitudeandBfieldDiff.Btot(:,1))+length(Modeldata.KvotAmplitudeandBfieldDiff.Btot(:,1)),:)=[ssave Modeldata.KvotAmplitudeandBfieldDiff.Btot];
    
   seig=s.*ones(length(Modeldata.Eigenvalues(:,1)),1); 
   Model.Eigenvalues(1+length(Model.Eigenvalues(:,1)):length(Model.Eigenvalues(:,1))+length(Modeldata.Eigenvalues(:,1)),:)=[seig Modeldata.Eigenvalues];
    
   %sB=s.*ones(length(Modeldata.Bfield.B1.Bx(:,1)),1); 
    %Model.Bfield.B1.Bx(1+length(Model.Bfield.B1.Bx(:,1)):length(Model.Bfield.B1.Bx(:,1))+length(Modeldata.Bfield.B1.Bx(:,1)),:)=[sB Modeldata.Bfield.B1.Bx];
   % Model.Bfield.B1.By(1+length(ModelModel.Bfield.B1.By(:,1)):length(Model.Bfield.B1.By(:,1))+length(Modeldata.Bfield.B1.By(:,1)),:)=[sB Modeldata.Bfield.B1.By];
   % Model.Bfield.B1.Bz(1+length(Model.Bfield.B1.Bz(:,1)):length(Model.Bfield.B1.Bz(:,1))+length(Modeldata.Bfield.B1.Bz(:,1)),:)=[sB Modeldata.Bfield.B1.Bz];
    %Model.Bfield.B2.Bx(1+length(Model.Bfield.B2.Bx(:,1)):length(Model.Bfield.B2.Bx(:,1))+length(Modeldata.Bfield.B2.Bx(:,1)),:)=[sB Modeldata.Bfield.B2.Bx];
    %Model.Bfield.B2.By(1+length(Model.Bfield.B2.By(:,1)):length(Model.Bfield.B2.By(:,1))+length(Modeldata.Bfield.B2.By(:,1)),:)=[sB Modeldata.Bfield.B2.By];
    %Model.Bfield.B2.Bz(1+length(Model.Bfield.B2.Bz(:,1)):length(Model.Bfield.B2.Bz(:,1))+length(Modeldata.Bfield.B2.Bz(:,1)),:)=[sB Modeldata.Bfield.B2.Bz];
    %Model.Bfield.B3.Bx(1+length(Model.Bfield.B3.Bx(:,1)):length(Model.Bfield.B3.Bx(:,1))+length(Modeldata.Bfield.B3.Bx(:,1)),:)=[sB Modeldata.Bfield.B3.Bx];
    %Model.Bfield.B3.By(1+length(Model.Bfield.B3.By(:,1)):length(Model.Bfield.B3.By(:,1))+length(Modeldata.Bfield.B3.By(:,1)),:)=[sB Modeldata.Bfield.B3.By];
    %Model.Bfield.B3.Bz(1+length(Model.Bfield.B3.Bz(:,1)):length(Model.Bfield.B3.Bz(:,1))+length(Modeldata.Bfield.B3.Bz(:,1)),:)=[sB Modeldata.Bfield.B3.Bz];
    %Model.Bfield.B4.Bx(1+length(Model.Bfield.B4.Bx(:,1)):length(Model.Bfield.B4.Bx(:,1))+length(Modeldata.Bfield.B4.Bx(:,1)),:)=[sB Modeldata.Bfield.B4.Bx];
    %Model.Bfield.B4.By(1+length(Model.Bfield.B4.By(:,1)):length(Model.Bfield.B4.By(:,1))+length(Modeldata.Bfield.B4.By(:,1)),:)=[sB Modeldata.Bfield.B4.By];
    %Model.Bfield.B4.Bz(1+length(Model.Bfield.B4.Bz(:,1)):length(Model.Bfield.B4.Bz(:,1))+length(Modeldata.Bfield.B4.Bz(:,1)),:)=[sB Modeldata.Bfield.B4.Bz];
    
    stype=s.*ones(length(Modeldata.NumberofType.A(:,1)),1); 
    Model.NumberofType.A(1+length(Model.NumberofType.A(:,1)):length(Model.NumberofType.A(:,1))+length(Modeldata.NumberofType.A(:,1)),:)=[stype Modeldata.NumberofType.A];
    Model.NumberofType.As(1+length(Model.NumberofType.As(:,1)):length(Model.NumberofType.As(:,1))+length(Modeldata.NumberofType.As(:,1)),:)=[stype Modeldata.NumberofType.As];
    Model.NumberofType.B(1+length(Model.NumberofType.B(:,1)):length(Model.NumberofType.B(:,1))+length(Modeldata.NumberofType.B(:,1)),:)=[stype Modeldata.NumberofType.B];
    Model.NumberofType.Bs(1+length(Model.NumberofType.Bs(:,1)):length(Model.NumberofType.Bs(:,1))+length(Modeldata.NumberofType.Bs(:,1)),:)=[stype Modeldata.NumberofType.Bs];
    Model.NumberofType.x(1+length(Model.NumberofType.x(:,1)):length(Model.NumberofType.x(:,1))+length(Modeldata.NumberofType.x(:,1)),:)=[stype Modeldata.NumberofType.x];
    Model.NumberofType.o(1+length(Model.NumberofType.o(:,1)):length(Model.NumberofType.o(:,1))+length(Modeldata.NumberofType.o(:,1)),:)=[stype Modeldata.NumberofType.o];
    Model.NumberofType.unknown(1+length(Model.NumberofType.unknown(:,1)):length(Model.NumberofType.unknown(:,1))+length(Modeldata.NumberofType.unknown(:,1)),:)=[stype Modeldata.NumberofType.unknown];
    
    serror=s.*ones(length(Modeldata.NumberofNullsError(:,1)),1); 
    Model.NumberofNullsError(1+length(Model.NumberofNullsError(:,1)):length(Model.NumberofNullsError(:,1))+length(Modeldata.NumberofNullsError(:,1)),:)=[serror Modeldata.NumberofNullsError];
    
    %sposerrorx=s.*ones(length(Modeldata.PositionError.x(:,1)),1);
    %sposerrory=s.*ones(length(Modeldata.PositionError.y(:,1)),1);
    %sposerrorz=s.*ones(length(Modeldata.PositionError.z(:,1)),1);
    
    %Model.PositionError.x(1+length(Model.PositionError.x(:,1)):length(Model.PositionError.x(:,1))+length(Modeldata.PositionError.x(:,1)),:)=[sposerrorx Modeldata.PositionError.x];
    %Model.PositionError.y(1+length(Model.PositionError.y(:,1)):length(Model.PositionError.y(:,1))+length(Modeldata.PositionError.y(:,1)),:)=[sposerrory Modeldata.PositionError.y];
    %Model.PositionError.z(1+length(Model.PositionError.z(:,1)):length(Model.PositionError.z(:,1))+length(Modeldata.PositionError.z(:,1)),:)=[sposerrorz Modeldata.PositionError.z];
end
%Increase (or decrease) in s for each ii run
s=s+50; %B %s=500;j=3;
%s=s-50; %A
end 
 disp('Saves data into datafile');
 save('ModelldataA.mat','Model','-v7.3');  
end
