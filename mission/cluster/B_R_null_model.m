function [Modeldata]=B_R_null_model(q,p,jperp,jparall,onesign,s,name)
%B_R_NULL_MODEL - Creates a simple model of linear B-field and position (R)
%taken from real data and calculates the number of different nulltypes
%found when a disturbance is added to a chosen satellite.
%p,q,s,jperp,parall and onesign are parameters from real data ...
%that has been rotated into the nulls coordinate system. See
%ROTATION_OF_GRADBDATA_TO_MODEL_PARAMETERS for more details

% jperp is the current perpendicular to the spine of the null and jparall is the
% current parallel to the spine of the null. s is the scaling parameter
% that is used to make the parameters unitless. onesign is the sign of the
% gradB(1,1) value. name is the name (in string) that you want to save the
% file as.

%Calculates the s/c's position
% % S/C 1-4 space coordinates using an estimated spacing between satellites
% from Cluster in Oct 1 2001
%Keeping the s/c's stuck and moving the null instead for simplicity sake
if nargin == 0
    help B_R_model;
    return;
elseif nargin < 7
    error('Too few input values. See usage: help B_R_null_model')
elseif nargin > 7
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

   jthresh=sqrt(((s*p)-(s*onesign)).^2+(s*q).^2);
   landa1=(s*p+s*onesign+sqrt((jthresh.^2-(s*jparall).^2)))/(2);
   landa2=(s*p+s*onesign-sqrt((jthresh.^2-(s*jparall).^2)))/(2);
   landa3=-1*(s*p+s*onesign);
   
   %Checks with the user which satellite should be disturbed
prompt = 'Which satellite do you want to disturb? 1/2/3/4: (default is 3)   ';
Sat = input(prompt,'s');
if strcmp(Sat,'')
    Sattest='3';
elseif strcmp(Sat,'1')
    Sattest=Sat;
elseif strcmp(Sat,'2')
    Sattest=Sat;
elseif strcmp(Sat,'3')
    Sattest=Sat;
elseif strcmp(Sat,'4')
    Sattest=Sat;    
else
    display('You did not a number between 1 and 4. Please restart the program and retry')
    return
end

%iii for loop runs to checks the distribution of the random function
for iii=1:250
    %Starting value for the amplitude of the disturbance
    amplitude=0;
    %i for loop loops over each amplitude increase so what i and iii does
    %toghether is calculate for each amplitude the distribution of the
    %random function since it will be changed for each iii run.
    for i=1:1000
        %The spherical coordinate angles
        phi=rand(size(tt)).*180;
        theta=rand(size(tt)).*360;
        
        
        %Disturbance in B-field of the chosen satellite
        Bz = amplitude.*cosd(phi);
        By = amplitude.*sind(theta).*sind(phi);
        Bx = amplitude.*cosd(theta).*sind(phi);
        
        
        
        %Magnetic field in nT measured at the first satellite
        B1x=s.*(onesign).*(x1-xn)+s.*((1/2).*(q-jparall).*(y1-yn));
        B1y=s.*(p).*(y1-yn)+s.*((1/2).*(q+jparall).*(x1-xn));
        B1z=s.*jperp.*(y1-yn)-s.*((p+onesign)).*(z1-zn);
        
        %Magnetic field in nT measured at the second satellite
        B2x=s.*(onesign).*(x2-xn)+s.*((1/2).*(q-jparall).*(y2-yn));
        B2y=s.*(p).*(y2-yn)+s.*((1/2).*(q+jparall).*(x2-xn));
        B2z=s.*jperp.*(y2-yn)-s.*((p+onesign)).*(z2-zn);
        
        %Magnetic field in nT measured at the third satellite
        B3x=s.*(onesign).*(x3-xn)+s.*((1/2).*(q-jparall).*(y3-yn));
        B3y=s.*(p).*(y3-yn)+s.*((1/2).*(q+jparall).*(x3-xn));
        B3z=s.*jperp.*(y3-yn)-s.*((p+onesign)).*(z3-zn);
        
        %Magnetic field in nT measured at the fourth satellite
        B4x=s.*(onesign).*(x4-xn)+s.*((1/2).*(q-jparall).*(y4-yn));
        B4y=s.*(p).*(y4-yn)+s.*((1/2).*(q+jparall).*(x4-xn));
        B4z=s.*jperp.*(y4-yn)-s.*((p+onesign)).*(z4-zn);
        
        if strcmp(Sattest,'1')
            B1x=B1x+Bx;
            B1y=B1y+By;
            B1z=B1z+Bz;
        end
        if strcmp(Sattest,'2')
            B2x=B2x+Bx;
            B2y=B2y+By;
            B2z=B2z+Bz;
        end
        if strcmp(Sattest,'3')
            B3x=B3x+Bx;
            B3y=B3y+By;
            B3z=B3z+Bz;
        end
        if strcmp(Sattest,'4')
            B4x=B4x+Bx;
            B4y=B4y+By;
            B4z=B4z+Bz;
        end
        
        
        B1=[tt B1x B1y B1z];
        B2=[tt B2x B2y B2z];
        B3=[tt B3x B3y B3z];
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
        nullposition=Nulls.nullPosition(nullPositionlogical,:);
        %saves only the Rn where we've located a null.
        %Rn(~nullPositionlogical,:)=NaN;
        
        %sigmaerror=amplitude.*ones(length(Nulls.nullPosition(:,2)),1);
        %xerror=[sigmaerror tt Rn(:,2)-Nulls.nullPosition(:,2)];
        %yerror=[sigmaerror tt Rn(:,3)-Nulls.nullPosition(:,3)];
        %zerror=[sigmaerror tt Rn(:,4)-Nulls.nullPosition(:,4)];
        
        %amplitudesave=amplitude.*ones(length(B1(:,2)),1);
        %time=B1(:,1);
        %Bfield1=[time amplitudesave B1(:,2:4)];
        %Bfield2=[time amplitudesave B2(:,2:4)];
        %Bfield3=[time amplitudesave B3(:,2:4)];
        %Bfield4=[time amplitudesave B4(:,2:4)];
        
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
        
        scalingamp=amplitude.*ones(length(nullposition(:,2)),1);
        GradB.Atype=[scalingamp Nulls.NullType.gradB.A(nullPositionlogical,2:end)];
        GradB.Astype=[scalingamp Nulls.NullType.gradB.As(nullPositionlogical,2:end)];
        GradB.Btype=[scalingamp Nulls.NullType.gradB.B(nullPositionlogical,2:end)];
        GradB.Bstype=[scalingamp Nulls.NullType.gradB.Bs(nullPositionlogical,2:end)];
        GradB.xtype=[scalingamp Nulls.NullType.gradB.x(nullPositionlogical,2:end)];
        GradB.otype=[scalingamp Nulls.NullType.gradB.o(nullPositionlogical,2:end)];
        GradB.unknowntype=[scalingamp Nulls.NullType.gradB.unknown(nullPositionlogical,2:end)];
        %Increases the amplitude with each i step so that in the end you will
        %have an amplitude which is about 75% of mBx. This is so the amplitude
        %is big enought for each increase of decrease when changing s in iii.
        amplitude=amplitude+((mBx(1,1)*0.75)/250);
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
            
            gradB.Atype=GradB.Atype;
            gradB.Btype=GradB.Btype;
            gradB.Astype=GradB.Astype;
            gradB.Bstype=GradB.Bstype;
            gradB.xtype=GradB.xtype;
            gradB.otype=GradB.otype;
            gradB.unknowntype=GradB.unknowntype;
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
           gradB.Atype(length(gradB.Atype(:,1))+1:length(GradB.Atype(:,1))+length(gradB.Atype(:,1)),:)=GradB.Atype;
        gradB.Btype(length(gradB.Btype(:,1))+1:length(GradB.Btype(:,1))+length(gradB.Btype(:,1)),:)=GradB.Btype;
        gradB.Astype(length( gradB.Astype(:,1))+1:length(GradB.Astype(:,1))+length(gradB.Astype(:,1)),:)=GradB.Astype;
        gradB.Bstype(length( gradB.Bstype(:,1))+1:length(GradB.Bstype(:,1))+length( gradB.Bstype(:,1)),:)=GradB.Bstype;
       gradB.xtype(length( gradB.xtype(:,1))+1:length(GradB.xtype(:,1))+length( gradB.xtype(:,1)),:)=GradB.xtype;
        gradB.otype(length( gradB.otype(:,1))+1:length(GradB.otype(:,1))+length( gradB.otype(:,1)),:)=GradB.otype;
        gradB.unknowntype(length( gradB.unknowntype(:,1))+1:length(GradB.unknowntype(:,1))+length( gradB.unknowntype(:,1)),:)=GradB.unknowntype;
            
            NumberofType.A(length(NumberofType.A(:,1))+1:length(A(:,1))+length(NumberofType.A(:,1)),:)=A;
            NumberofType.As(length(NumberofType.As(:,1))+1:length(As(:,1))+length(NumberofType.As(:,1)),:)=As;
            NumberofType.B(length(NumberofType.B(:,1))+1:length(B(:,1))+length(NumberofType.B(:,1)),:)=B;
            NumberofType.Bs(length(NumberofType.Bs(:,1))+1:length(Bs(:,1))+length(NumberofType.Bs(:,1)),:)=Bs;
            NumberofType.x(length(NumberofType.x(:,1))+1:length(x(:,1))+length(NumberofType.x(:,1)),:)=x;
            NumberofType.o(length(NumberofType.o(:,1))+1:length(o(:,1))+length(NumberofType.o(:,1)),:)=o;
            NumberofType.unknown(length(NumberofType.unknown(:,1))+1:length(unknown(:,1))+length(NumberofType.unknown(:,1)),:)=unknown;
            
            NumberofNullsError(length(NumberofNullsError(:,1))+1:length(numberofnullserror(:,1))+length(NumberofNullsError(:,1)),:)=numberofnullserror;
            
            %PositionError.x(length(PositionError.x(:,1))+1:length(xerror(:,1))+length(PositionError.x(:,1)),:)=xerror;
            %PositionError.y(length(PositionError.y(:,1))+1:length(yerror(:,1))+length(PositionError.y(:,1)),:)=yerror;
            %PositionError.z(length(PositionError.z(:,1))+1:length(zerror(:,1))+length(PositionError.z(:,1)),:)=zerror;
        end
    end
    display('Loop nr: ', num2str(iii))
    %Saving data for each iii run (random function distribution)
    if iii==1
        Modeldata.KvotAmplitudeandBfieldDiff.Bx=KvotAmplitudeandBfieldDiff.Bx;
        Modeldata.KvotAmplitudeandBfieldDiff.By=KvotAmplitudeandBfieldDiff.By;
        Modeldata.KvotAmplitudeandBfieldDiff.Bz=KvotAmplitudeandBfieldDiff.Bz;
        Modeldata.KvotAmplitudeandBfieldDiff.Btot=KvotAmplitudeandBfieldDiff.Btot;
        
         Modeldata.GradB.A=gradB.Atype;
         Modeldata.GradB.B=gradB.Btype;
        Modeldata.GradB.As=gradB.Astype;
        Modeldata.GradB.Bs=gradB.Bstype;
            Modeldata.GradB.x=gradB.xtype;
            Modeldata.GradB.o=gradB.otype;
           Modeldata.GradB.unknown=gradB.unknowntype;
        
        %Modeldata.Bfield.B1.Bx=[Bfield.B1(:,1) Bfield.B1(:,2) Bfield.B1(:,3)];
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
        
        Modeldata.scale.p=p;
        Modeldata.scale.q=q;
        Modeldata.scale.s=s;
        Modeldata.scale.jthresh=jthresh;
        Modeldata.scale.jparall=jparall;
        Modeldata.scale.jperp=jperp;
        Modeldata.scale.onesign=onesign;
        
        
        %Modeldata.PositionError.x=PositionError.x;
       % Modeldata.PositionError.y=PositionError.y;
        %Modeldata.PositionError.z=PositionError.z;
        
        %Saves data for each next iii step in the next column for random function
        %distribution
    else
        disp('Saves data into Modeldata file');
        Modeldata.KvotAmplitudeandBfieldDiff.Bx(:,length(Modeldata.KvotAmplitudeandBfieldDiff.Bx(1,:))+1:1+length(Modeldata.KvotAmplitudeandBfieldDiff.Bx(1,:)))=KvotAmplitudeandBfieldDiff.Bx(:,2);
        Modeldata.KvotAmplitudeandBfieldDiff.By(:,length(Modeldata.KvotAmplitudeandBfieldDiff.By(1,:))+1:1+length(Modeldata.KvotAmplitudeandBfieldDiff.By(1,:)))=KvotAmplitudeandBfieldDiff.By(:,2);
        Modeldata.KvotAmplitudeandBfieldDiff.Bz(:,length(Modeldata.KvotAmplitudeandBfieldDiff.Bz(1,:))+1:1+length(Modeldata.KvotAmplitudeandBfieldDiff.Bz(1,:)))=KvotAmplitudeandBfieldDiff.Bz(:,2);
        Modeldata.KvotAmplitudeandBfieldDiff.Btot(:,length(Modeldata.KvotAmplitudeandBfieldDiff.Btot(1,:))+1:1+length(Modeldata.KvotAmplitudeandBfieldDiff.Btot(1,:)))=KvotAmplitudeandBfieldDiff.Btot(:,2);
        
            
        %Modeldata.Bfield.B1.Bx(:,length(Modeldata.Bfield.B1.Bx(1,:))+1:1+length(Modeldata.Bfield.B1.Bx(1,:)))=Bfield.B1(:,3);
        %Modeldata.Bfield.B1.By(:,length(Modeldata.Bfield.B1.By(1,:))+1:1+length(Modeldata.Bfield.B1.By(1,:)))=Bfield.B1(:,4);
        %Modeldata.Bfield.B1.Bz(:,length(Modeldata.Bfield.B1.Bz(1,:))+1:1+length(Modeldata.Bfield.B1.Bz(1,:)))=Bfield.B1(:,5);
        %Modeldata.Bfield.B2.Bx(:,length(Modeldata.Bfield.B2.Bx(1,:))+1:1+length(Modeldata.Bfield.B2.Bx(1,:)))=Bfield.B2(:,3);
        %Modeldata.Bfield.B2.By(:,length(Modeldata.Bfield.B2.By(1,:))+1:1+length(Modeldata.Bfield.B2.By(1,:)))=Bfield.B2(:,4);
        %Modeldata.Bfield.B2.Bz(:,length(Modeldata.Bfield.B2.Bz(1,:))+1:1+length(Modeldata.Bfield.B2.Bz(1,:)))=Bfield.B2(:,5);
        %Modeldata.Bfield.B3.Bx(:,length(Modeldata.Bfield.B3.Bx(1,:))+1:1+length(Modeldata.Bfield.B3.Bx(1,:)))=Bfield.B3(:,3);
        %Modeldata.Bfield.B3.By(:,length(Modeldata.Bfield.B3.By(1,:))+1:1+length(Modeldata.Bfield.B3.By(1,:)))=Bfield.B3(:,4);
        %Modeldata.Bfield.B3.Bz(:,length(Modeldata.Bfield.B3.Bz(1,:))+1:1+length(Modeldata.Bfield.B3.Bz(1,:)))=Bfield.B3(:,5);
        %Modeldata.Bfield.B4.Bx(:,length(Modeldata.Bfield.B4.Bx(1,:))+1:1+length(Modeldata.Bfield.B4.Bx(1,:)))=Bfield.B4(:,3);
        %Modeldata.Bfield.B4.By(:,length(Modeldata.Bfield.B4.By(1,:))+1:1+length(Modeldata.Bfield.B4.By(1,:)))=Bfield.B4(:,4);
        %Modeldata.Bfield.B4.Bz(:,length(Modeldata.Bfield.B4.Bz(1,:))+1:1+length(Modeldata.Bfield.B4.Bz(1,:)))=Bfield.B4(:,5);
        
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
end
%Saves data for each run with different s values.
disp('Saves data into datafile');
save([name,'.mat'],'Modeldata','-v7.3');
end
