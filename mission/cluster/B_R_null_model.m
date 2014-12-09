function [Modeldata]=B_R_null_model(q,p,jperp,jparall)
%B_R_null_model - Creates a simple model of B-field with no currents and position (R)
%measurements for 4 s/c's.

% This function creates a simple linear Bfield that can have either complex
% or real eigenvalues. p,q are parameters from real data.
% jperp is the current perpendicular to the spine and jparall is the
% current parallel to the spine.

% Important: Time tags should be included with the vectors
%
%

%Calculates the s/c's position
% % S/C 1-4 space coordinates using an estimated spacing between satellites
% from Cluster in Oct 1 2001
%Keeping the s/c's stuck and moving the null instead for simplicity sake
if nargin == 0
    help B_R_model;
    return;
elseif nargin < 4
    error('Too few input values. See usage: help B_R_model')
elseif nargin > 4
    error('Too many input values. See usage: help B_R_model')
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
    x1(i,1)=-106493;
    x2(i,1)=-106558.601562500;
    x3(i,1)=-106674.898437500;
    x4(i,1)=-106561.571555391;
    y1(i,1)=-58033.9495605081;
    y2(i,1)=-57835.8976059652;
    y3(i,1)=-57916.6523337731;
    y4(i,1)=-58032.3476151215;
    z1(i,1)=2951.23366036776;
    z2(i,1)=2855.83738888417;
    z3(i,1)=2929.01943807376;
    z4(i,1)=2762.29629071977;
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
x0=-106700;
y0=-58100;
z0=2700;
xn=x0+20.*t;
yn=y0+20.*t;
zn=z0+20.*t;
Rn=[tt xn yn zn];
%Checks with the user if the inserted data is in the solar wind
prompt = 'Do you want a A(As) type null? yes/no: (default is yes)   ';
Type = input(prompt,'s');
if strcmp(Type,'')
    typetest='yes';
elseif strcmp(Type,'no') || strcmp(Type,'yes')
    typetest=Type;
else
    display('You did not pick yes/no. Please restart the program')
    return
end
%Threshold current
%Eigenvalues
if strcmp(typetest,'yes')
    jthresh=sqrt((p+1).^2+q.^2);
landa1=(p-1+sqrt(jthresh.^2-jparall.^2))/(2);
landa2=(p-1-sqrt(jthresh.^2-jparall.^2))/(2);
landa3=-1*(p-1);
else
    jthresh=sqrt((p-1).^2+q.^2);
landa1=(p+1+sqrt(jthresh.^2-jparall.^2))/(2);
landa2=(p+1-sqrt(jthresh.^2-jparall.^2))/(2);
   landa3=-1*(p+1); 
end
%iii for loop runs to checks the distribution of the random function
for iii=1:500
    %Starting value for the amplitude of the disturbance
    amplitude=0;
    uncertaincyamplitude1=0;
    uncertaincyamplitude2=0;
    uncertaincyamplitude4=0;
    %i for loop loops over each amplitude increase so what i and iii does
    %toghether is calculate for each amplitude the distribution of the
    %random function since it will be changed for each iii run.
    for i=1:1000
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
        
        if strcmp(typetest,'yes')
   
        %Magnetic field measured at the third satellite
        B3x=(-1).*(x3-xn)+((1/2).*(q+jparall).*(y3-yn));
        B3y=(p).*(y3-yn)+((1/2).*(q-jparall).*(x3-xn));
        B3z=jperp.*(y3-yn)-((p-1)).*(z3-zn);
        B3x=B3x+Bx;
        B3y=B3y+By;
        B3z=B3z+Bz;
        B3=[tt B3x B3y B3z];
        
        %Magnetic field measured at the first satellite
        B1x=(-1).*(x1-xn)+((1/2).*(q+jparall).*(y1-yn));
        B1y=(p).*(y1-yn)+((1/2).*(q-jparall).*(x1-xn));
        B1z=jperp.*(y1-yn)-((p-1)).*(z1-zn);
        B1x=B1x+Buncertaintyx1;
        B1y=B1y+Buncertaintyy1;
        B1z=B1z+Buncertaintyz1;
        B1=[tt B1x B1y B1z];
        
        %Magnetic field measured at the second satellite
        B2x=(-1).*(x2-xn)+((1/2).*(q+jparall).*(y2-yn));
        B2y=(p).*(y2-yn)+((1/2).*(q-jparall).*(x2-xn));
        B2z=jperp.*(y2-yn)-((p-1)).*(z2-zn);
        B2x=B2x+Buncertaintyx2;
        B2y=B2y+Buncertaintyy2;
        B2z=B2z+Buncertaintyz2;
        B2=[tt B2x B2y B2z];
        
        %Magnetic field measured at the fourth satellite
        B4x=(-1).*(x4-xn)+((1/2).*(q+jparall).*(y4-yn));
        B4y=(p).*(y4-yn)+((1/2).*(q-jparall).*(x4-xn));
        B4z=jperp.*(y4-yn)-((p-1)).*(z4-zn);
        B4x=B4x+Buncertaintyx4;
        B4y=B4y+Buncertaintyy4;
        B4z=B4z+Buncertaintyz4;
        B4=[tt B4x B4y B4z];
        
        else
             %Magnetic field measured at the third satellite
        B3x=(1).*(x3-xn)+((1/2).*(q-jparall).*(y3-yn));
        B3y=(p).*(y3-yn)+((1/2).*(q+jparall).*(x3-xn));
        B3z=jperp.*(y3-yn)-((p+1)).*(z3-zn);
        B3x=B3x+Bx;
        B3y=B3y+By;
        B3z=B3z+Bz;
        B3=[tt B3x B3y B3z];
        
        %Magnetic field measured at the first satellite
        B1x=(1).*(x1-xn)+((1/2).*(q-jparall).*(y1-yn));
        B1y=(p).*(y1-yn)+((1/2).*(q+jparall).*(x1-xn));
        B1z=jperp.*(y1-yn)-((p+1)).*(z1-zn);
        B1x=B1x+Buncertaintyx1;
        B1y=B1y+Buncertaintyy1;
        B1z=B1z+Buncertaintyz1;
        B1=[tt B1x B1y B1z];
        
        %Magnetic field measured at the second satellite
        B2x=(1).*(x2-xn)+((1/2).*(q-jparall).*(y2-yn));
        B2y=(p).*(y2-yn)+((1/2).*(q+jparall).*(x2-xn));
        B2z=jperp.*(y2-yn)-((p+1)).*(z2-zn);
        B2x=B2x+Buncertaintyx2;
        B2y=B2y+Buncertaintyy2;
        B2z=B2z+Buncertaintyz2;
        B2=[tt B2x B2y B2z];
        
        %Magnetic field measured at the fourth satellite
        B4x=(1).*(x4-xn)+((1/2).*(q-jparall).*(y4-yn));
        B4y=(p).*(y4-yn)+((1/2).*(q+jparall).*(x4-xn));
        B4z=jperp.*(y4-yn)-((p+1)).*(z4-zn);
        B4x=B4x+Buncertaintyx4;
        B4y=B4y+Buncertaintyy4;
        B4z=B4z+Buncertaintyz4;
        B4=[tt B4x B4y B4z];
        end
        
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
        Nulls=c_4_null(R1,R2,R3,R4,B1,B2,B3,B4);
        
        %Here starts the saving function
        nullPositionlogical = false(length(Nulls.nullPosition(:,1)),1);
    for TimeInterval2=1:length(Nulls.nullPosition(:,1))
        if all(isnan(Nulls.nullPosition(TimeInterval2,2:4))) %If the nullPosition has been NaN continue to the next nullPosition
            continue
        end
        nullPositionlogical(TimeInterval2,1)=true;
    end
        disp('calculates number of type etc.');
        % nullPosition(~nullPositionlogical,:)=NaN;
        %saves only the Rn where we've located a null.
        Rn(~nullPositionlogical,:)=NaN;
        
        sigmaerror=amplitude.*ones(length(Nulls.nullPosition(:,2)),1);
        xerror=[sigmaerror tt Rn(:,2)-Nulls.nullPosition(:,2)];
        yerror=[sigmaerror tt Rn(:,3)-Nulls.nullPosition(:,3)];
        zerror=[sigmaerror tt Rn(:,4)-Nulls.nullPosition(:,4)];
        
        amplitudesave=amplitude.*ones(length(B1(:,2)),1);
        time=B1(:,1);
        Bfield1=[time amplitudesave B1(:,2:4)];
        Bfield2=[time amplitudesave B2(:,2:4)];
        Bfield3=[time amplitudesave B3(:,2:4)];
        Bfield4=[time amplitudesave B4(:,2:4)];
        
        %Saves the number of each type found for each amplitude 500 is the size
        %of the Bfields so I remove all NaN values where a certain type hasn't
        %been found.
        A=[amplitude 500-sum(isnan(Nulls.NullType.position.A(:,2)))];
        B=[amplitude 500-sum(isnan(Nulls.NullType.position.B(:,2)))];
        As=[amplitude 500-sum(isnan(Nulls.NullType.position.As(:,2)))];
        Bs=[amplitude 500-sum(isnan(Nulls.NullType.position.Bs(:,2)))];
        unknown=[amplitude 500-sum(isnan(Nulls.NullType.position.unknown(:,2)))];
        x=[amplitude 500-sum(isnan(Nulls.NullType.position.x(:,2)))];
        o=[amplitude 500-sum(isnan(Nulls.NullType.position.o(:,2)))];
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
            
            Bfield.B1=Bfield1;
            Bfield.B2=Bfield2;
            Bfield.B3=Bfield3;
            Bfield.B4=Bfield4;
            
            NumberofType.A=A;
            NumberofType.As=As;
            NumberofType.B=B;
            NumberofType.Bs=Bs;
            NumberofType.unknown=unknown;
            NumberofType.x=x;
            NumberofType.o=o;
            
            NumberofNullsError=numberofnullserror;
            
            PositionError.x=xerror;
            PositionError.y=yerror;
            PositionError.z=zerror;
            
            %Saves data for each next i step in the next row for the different
            %amplitudes
        else
            KvotAmplitudeandBfieldDiff.Bx(length(KvotAmplitudeandBfieldDiff.Bx(:,1))+1:length(KvotBetweenamplitudeandBxieldDifference(:,1))+length(KvotAmplitudeandBfieldDiff.Bx(:,1)),:)=KvotBetweenamplitudeandBxieldDifference;
            KvotAmplitudeandBfieldDiff.By(length(KvotAmplitudeandBfieldDiff.By(:,1))+1:length(KvotBetweenamplitudeandByieldDifference(:,1))+length(KvotAmplitudeandBfieldDiff.By(:,1)),:)=KvotBetweenamplitudeandByieldDifference;
            KvotAmplitudeandBfieldDiff.Bz(length(KvotAmplitudeandBfieldDiff.Bz(:,1))+1:length(KvotBetweenamplitudeandBzieldDifference(:,1))+length(KvotAmplitudeandBfieldDiff.Bz(:,1)),:)=KvotBetweenamplitudeandBzieldDifference;
            KvotAmplitudeandBfieldDiff.Btot(length(KvotAmplitudeandBfieldDiff.Btot(:,1))+1:length(KvotBetweenamplitudeandBtotDifference(:,1))+length(KvotAmplitudeandBfieldDiff.Btot(:,1)),:)=KvotBetweenamplitudeandBtotDifference;
            
            Bfield.B1(length(Bfield.B1(:,1))+1:length(Bfield1(:,1))+length(Bfield.B1(:,1)),:)=Bfield1;
            Bfield.B2(length(Bfield.B2(:,1))+1:length(Bfield2(:,1))+length(Bfield.B2(:,1)),:)=Bfield2;
            Bfield.B3(length(Bfield.B3(:,1))+1:length(Bfield3(:,1))+length(Bfield.B3(:,1)),:)=Bfield3;
            Bfield.B4(length(Bfield.B4(:,1))+1:length(Bfield4(:,1))+length(Bfield.B4(:,1)),:)=Bfield4;
            
            NumberofType.A(length(NumberofType.A(:,1))+1:length(A(:,1))+length(NumberofType.A(:,1)),:)=A;
            NumberofType.As(length(NumberofType.As(:,1))+1:length(As(:,1))+length(NumberofType.As(:,1)),:)=As;
            NumberofType.B(length(NumberofType.B(:,1))+1:length(B(:,1))+length(NumberofType.B(:,1)),:)=B;
            NumberofType.Bs(length(NumberofType.Bs(:,1))+1:length(Bs(:,1))+length(NumberofType.Bs(:,1)),:)=Bs;
            NumberofType.x(length(NumberofType.x(:,1))+1:length(x(:,1))+length(NumberofType.x(:,1)),:)=x;
            NumberofType.o(length(NumberofType.o(:,1))+1:length(o(:,1))+length(NumberofType.o(:,1)),:)=o;
            NumberofType.unknown(length(NumberofType.unknown(:,1))+1:length(unknown(:,1))+length(NumberofType.unknown(:,1)),:)=unknown;
            
            NumberofNullsError(length(NumberofNullsError(:,1))+1:length(numberofnullserror(:,1))+length(NumberofNullsError(:,1)),:)=numberofnullserror;
            
            PositionError.x(length(PositionError.x(:,1))+1:length(xerror(:,1))+length(PositionError.x(:,1)),:)=xerror;
            PositionError.y(length(PositionError.y(:,1))+1:length(yerror(:,1))+length(PositionError.y(:,1)),:)=yerror;
            PositionError.z(length(PositionError.z(:,1))+1:length(zerror(:,1))+length(PositionError.z(:,1)),:)=zerror;
        end
    end
    display('Loop nr: ', num2str(iii))
    %Saving data for each iii run (random function distribution)
    if iii==1
        Modeldata.KvotAmplitudeandBfieldDiff.Bx=KvotAmplitudeandBfieldDiff.Bx;
        Modeldata.KvotAmplitudeandBfieldDiff.By=KvotAmplitudeandBfieldDiff.By;
        Modeldata.KvotAmplitudeandBfieldDiff.Bz=KvotAmplitudeandBfieldDiff.Bz;
        Modeldata.KvotAmplitudeandBfieldDiff.Btot=KvotAmplitudeandBfieldDiff.Btot;
        
        Modeldata.Bfield.B1.Bx=[Bfield.B1(:,1) Bfield.B1(:,2) Bfield.B1(:,3)];
        Modeldata.Bfield.B1.By=[Bfield.B1(:,1) Bfield.B1(:,2) Bfield.B1(:,4)];
        Modeldata.Bfield.B1.Bz=[Bfield.B1(:,1) Bfield.B1(:,2) Bfield.B1(:,5)];
        Modeldata.Bfield.B2.Bx=[Bfield.B2(:,1) Bfield.B2(:,2) Bfield.B2(:,3)];
        Modeldata.Bfield.B2.By=[Bfield.B2(:,1) Bfield.B2(:,2) Bfield.B2(:,4)];
        Modeldata.Bfield.B2.Bz=[Bfield.B2(:,1) Bfield.B2(:,2) Bfield.B2(:,5)];
        Modeldata.Bfield.B3.Bx=[Bfield.B3(:,1) Bfield.B3(:,2) Bfield.B3(:,3)];
        Modeldata.Bfield.B3.By=[Bfield.B3(:,1) Bfield.B3(:,2) Bfield.B3(:,4)];
        Modeldata.Bfield.B3.Bz=[Bfield.B3(:,1) Bfield.B3(:,2) Bfield.B3(:,5)];
        Modeldata.Bfield.B4.Bx=[Bfield.B4(:,1) Bfield.B4(:,2) Bfield.B4(:,3)];
        Modeldata.Bfield.B4.By=[Bfield.B4(:,1) Bfield.B4(:,2) Bfield.B4(:,4)];
        Modeldata.Bfield.B4.Bz=[Bfield.B4(:,1) Bfield.B4(:,2) Bfield.B4(:,5)];
        
        Modeldata.NumberofType.A=NumberofType.A;
        Modeldata.NumberofType.As=NumberofType.As;
        Modeldata.NumberofType.B=NumberofType.B;
        Modeldata.NumberofType.Bs=NumberofType.Bs;
        Modeldata.NumberofType.unknown=NumberofType.unknown;
        Modeldata.NumberofType.x=NumberofType.x;
        Modeldata.NumberofType.o=NumberofType.o;
        
        Modeldata.Eigenvalues=[landa1 landa2 landa3];
        
        Modeldata.NumberofNullsError=NumberofNullsError;
        
        Modeldata.scale.p=p;
        Modeldata.scale.q=q;
        Modeldata.scale.jparall=jparall;
        Modeldata.scale.jperp=jperp;
        
        
        Modeldata.PositionError.x=PositionError.x;
        Modeldata.PositionError.y=PositionError.y;
        Modeldata.PositionError.z=PositionError.z;
        
        %Saves data for each next iii step in the next column for random function
        %distribution
    else
        disp('Saves data into Modeldata file');
        Modeldata.KvotAmplitudeandBfieldDiff.Bx(:,length(Modeldata.KvotAmplitudeandBfieldDiff.Bx(1,:))+1:1+length(Modeldata.KvotAmplitudeandBfieldDiff.Bx(1,:)))=KvotAmplitudeandBfieldDiff.Bx(:,2);
        Modeldata.KvotAmplitudeandBfieldDiff.By(:,length(Modeldata.KvotAmplitudeandBfieldDiff.By(1,:))+1:1+length(Modeldata.KvotAmplitudeandBfieldDiff.By(1,:)))=KvotAmplitudeandBfieldDiff.By(:,2);
        Modeldata.KvotAmplitudeandBfieldDiff.Bz(:,length(Modeldata.KvotAmplitudeandBfieldDiff.Bz(1,:))+1:1+length(Modeldata.KvotAmplitudeandBfieldDiff.Bz(1,:)))=KvotAmplitudeandBfieldDiff.Bz(:,2);
        Modeldata.KvotAmplitudeandBfieldDiff.Btot(:,length(Modeldata.KvotAmplitudeandBfieldDiff.Btot(1,:))+1:1+length(Modeldata.KvotAmplitudeandBfieldDiff.Btot(1,:)))=KvotAmplitudeandBfieldDiff.Btot(:,2);
        
        Modeldata.Bfield.B1.Bx(:,length(Modeldata.Bfield.B1.Bx(1,:))+1:1+length(Modeldata.Bfield.B1.Bx(1,:)))=Bfield.B1(:,3);
        Modeldata.Bfield.B1.By(:,length(Modeldata.Bfield.B1.By(1,:))+1:1+length(Modeldata.Bfield.B1.By(1,:)))=Bfield.B1(:,4);
        Modeldata.Bfield.B1.Bz(:,length(Modeldata.Bfield.B1.Bz(1,:))+1:1+length(Modeldata.Bfield.B1.Bz(1,:)))=Bfield.B1(:,5);
        Modeldata.Bfield.B2.Bx(:,length(Modeldata.Bfield.B2.Bx(1,:))+1:1+length(Modeldata.Bfield.B2.Bx(1,:)))=Bfield.B2(:,3);
        Modeldata.Bfield.B2.By(:,length(Modeldata.Bfield.B2.By(1,:))+1:1+length(Modeldata.Bfield.B2.By(1,:)))=Bfield.B2(:,4);
        Modeldata.Bfield.B2.Bz(:,length(Modeldata.Bfield.B2.Bz(1,:))+1:1+length(Modeldata.Bfield.B2.Bz(1,:)))=Bfield.B2(:,5);
        Modeldata.Bfield.B3.Bx(:,length(Modeldata.Bfield.B3.Bx(1,:))+1:1+length(Modeldata.Bfield.B3.Bx(1,:)))=Bfield.B3(:,3);
        Modeldata.Bfield.B3.By(:,length(Modeldata.Bfield.B3.By(1,:))+1:1+length(Modeldata.Bfield.B3.By(1,:)))=Bfield.B3(:,4);
        Modeldata.Bfield.B3.Bz(:,length(Modeldata.Bfield.B3.Bz(1,:))+1:1+length(Modeldata.Bfield.B3.Bz(1,:)))=Bfield.B3(:,5);
        Modeldata.Bfield.B4.Bx(:,length(Modeldata.Bfield.B4.Bx(1,:))+1:1+length(Modeldata.Bfield.B4.Bx(1,:)))=Bfield.B4(:,3);
        Modeldata.Bfield.B4.By(:,length(Modeldata.Bfield.B4.By(1,:))+1:1+length(Modeldata.Bfield.B4.By(1,:)))=Bfield.B4(:,4);
        Modeldata.Bfield.B4.Bz(:,length(Modeldata.Bfield.B4.Bz(1,:))+1:1+length(Modeldata.Bfield.B4.Bz(1,:)))=Bfield.B4(:,5);
        
        Modeldata.NumberofType.A(:,length(Modeldata.NumberofType.A(1,:))+1:1+length(Modeldata.NumberofType.A(1,:)))=NumberofType.A(:,2);
        Modeldata.NumberofType.As(:,length(Modeldata.NumberofType.As(1,:))+1:1+length(Modeldata.NumberofType.As(1,:)))=NumberofType.As(:,2);
        Modeldata.NumberofType.B(:,length(Modeldata.NumberofType.B(1,:))+1:1+length(Modeldata.NumberofType.B(1,:)))=NumberofType.B(:,2);
        Modeldata.NumberofType.Bs(:,length(Modeldata.NumberofType.Bs(1,:))+1:1+length(Modeldata.NumberofType.Bs(1,:)))=NumberofType.Bs(:,2);
        Modeldata.NumberofType.x(:,length(Modeldata.NumberofType.x(1,:))+1:1+length(Modeldata.NumberofType.x(1,:)))=NumberofType.x(:,2);
        Modeldata.NumberofType.o(:,length(Modeldata.NumberofType.o(1,:))+1:1+length(Modeldata.NumberofType.o(1,:)))=NumberofType.o(:,2);
        Modeldata.NumberofType.unknown(:,length(Modeldata.NumberofType.unknown(1,:))+1:1+length(Modeldata.NumberofType.unknown(1,:)))=NumberofType.unknown(:,2);
        
        Modeldata.NumberofNullsError(:,length(Modeldata.NumberofNullsError(1,:))+1:1+length(Modeldata.NumberofNullsError(1,:)))=NumberofNullsError(:,2);
        
        Modeldata.PositionError.x(:,length(Modeldata.PositionError.x(1,:))+1:1+length(Modeldata.PositionError.x(1,:)))=PositionError.x(:,3);
        Modeldata.PositionError.y(:,length(Modeldata.PositionError.y(1,:))+1:1+length(Modeldata.PositionError.y(1,:)))=PositionError.y(:,3);
        Modeldata.PositionError.z(:,length(Modeldata.PositionError.z(1,:))+1:1+length(Modeldata.PositionError.z(1,:)))=PositionError.z(:,3);
    end
end

 disp('Saves data into datafile');
 name='Modeldata.mat';
save(name,'Modeldata','-v7.3'); 
end
