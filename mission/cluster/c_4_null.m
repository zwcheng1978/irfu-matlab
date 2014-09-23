function Nulls=c_4_null(R1,R2,R3,R4,B1,B2,B3,B4,varargin)
%C_4_NULL_POSITION - Calculates the null position within the tetrahedron using 4 spacecraft technique
%
%This function calculates the null position within the tetrahedron made up 
%of the satellites by assuming that the B-field in the vicinity of the 
%null can be determined by considering a Taylor
%expansion of the lowest order of B about the null.
%
%   [nullPosition,C4limits,dRmin,NullType,Requirement]=C_4_NULL_POSITION(R1,R2,R3,R4,B1,B2,B3,B4);
%   [nullPosition,C4limits,dRmin,NullType,Requirement]=C_4_NULL_POSITION(R1,R2,R3,R4,B1,B2,B3,B4, threshold_value,scseparation_value);
%   -threshold=100 means no restriction. Default value for the satellite
%   separation is 1000km. The value needs to be given in km. If a threshold
%   or sc separation value is given you need to also give the other value
%   for the program to work.
%   OUTPUT
%   nullPosition = [Time xn yn zn]
%   Requirement is a structure containing the two restrictions used on the
%   data. BfieldandEigenvalueslessthanchosenpercentage is if the
%   eigenvalues and Bfield values for the datapoints are below the chosen
%   percentage value (default=40%). DistancewithinSCconfiguration is the
%   logical vector showing true (1) if the null point is within the SC
%   configuration.
%   dRmin = [Time dRmin] Gives the minimum distance to the null point
%   looking from all satellites
%   C4limits is a structure containing the maxmimum and minimum positions
%   of all satellites at each time tag.
%   NullType is a structure containing two structures where each nulltype identified at each
%   time tag is given in the minimum distance from the satellites to the nullpoint and in the exact position (x,y,z).
%   Eigenvalues is a structure consisting of the lambda values and another
%   structure with the eigenvalues of the different lamda values. Lambda 1
%   of the lambda values coincide with the 3 first values on the first row.
%   Lambda 2 of the lambda values coincide with the 3 middle values on the first row.
% INPUT
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
if length(varargin)==1
    error('Too few input values. See usage: help c_4_null_position')
end
if isempty(varargin) == true
    threshold = 40;
    scseparation=1000; %Ion intertial length
else
    threshold = varargin{1};
    scseparation=varargin{2};
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
% Check if there's any nulls inside the s/c tetrahedron using Poincaré
% index


%Check which type the nulls are
[NullsEigen,Eigenvaluestypes,Current,constraint,errors]=c_4_null_type(R1,R2,R3,R4,B1,B2,B3,B4,threshold);

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
dRmin(:,1) = Time;

%Null position
Rn1 = irf_add(1,R1,-1,dR1);  %R1-dR1
Rn2 = irf_add(1,R2,-1,dR2);
Rn3 = irf_add(1,R3,-1,dR3);
Rn4 = irf_add(1,R4,-1,dR4);
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


%For each eigenvalue corresponding to the tolerance level (the two errors less or equal to 40%) break out their corresponding time and dR value
%(the minimum distance from all s/c to the null)
disp('Sorting based on the null located within the s/c tetrahedron thus having a maximum length of one ion inertial length');
sortNull=dRmin(:,2) <= scseparation;
sortNullDx=Rn1(:,2) >= minX & Rn1(:,2) <= maxX;
sortNullDy=Rn1(:,3) >= minY & Rn1(:,3) <= maxY;
sortNullDz=Rn1(:,4) >= minZ & Rn1(:,4) <= maxZ;
sortdr= sortNull & sortNullDx & sortNullDy & sortNullDz;

Nulls.errors=errors;
Nulls.index(~constraint,:)=NaN;
Nulls.index(~sortdr,:)=NaN;

Current.jParallel(~sortdr,:)=NaN;
Nulls.current.jParallel=Current.jParallel;

Current.jPerpendicular(~sortdr,:)=NaN;
Nulls.current.jPerpendicular=Current.jPerpendicular;

% min and max for all s/c's
minX(~constraint,:) = NaN;
maxX(~constraint,:) = NaN;
minY(~constraint,:) = NaN;
maxY(~constraint,:) = NaN;
minZ(~constraint,:) = NaN;
maxZ(~constraint,:) = NaN;

%Spacecraft max,min position
Nulls.C4limits.minX = minX;
Nulls.C4limits.maxX = maxX;
Nulls.C4limits.minY = minY;
Nulls.C4limits.maxY = maxY;
Nulls.C4limits.minZ = minZ;
Nulls.C4limits.maxZ = maxZ;

%The null positions
Time=B1(:,1);
Rn(~constraint,:) = NaN;
Rn(~sortdr,:)     = NaN;
Nulls.nullPosition      = [Time Rn];

Rn1(~constraint,:) = NaN;
Rn1(~sortdr,:)     = NaN;
Rn2(~constraint,:) = NaN;
Rn2(~sortdr,:)     = NaN;

Rn3(~constraint,:) = NaN;
Rn3(~sortdr,:)     = NaN;

Rn4(~constraint,:) = NaN;
Rn4(~sortdr,:)     = NaN;

%dRmin
dRmin(~constraint,2)= NaN;

%Eigenvalues
Eigenvaluestypes.lambda(~sortdr,:)     = NaN;
Eigenvaluestypes.eigenvectors.lambda1(~sortdr,:)     = NaN;
Eigenvaluestypes.eigenvectors.lambda2(~sortdr,:)     = NaN;
Eigenvaluestypes.eigenvectors.lambda3(~sortdr,:)     = NaN;
ttt                        =Time;
ttt(~constraint,:)         = NaN;
ttt(~sortdr,:)             = NaN;


EigenvaluesA               = Eigenvaluestypes.lambda;
EigenvaluesA(~NullsEigen.eigA,:)= NaN;
tA                         =ttt;
tA(~NullsEigen.eigA,1)          =NaN;
Nulls.Eigenvalues.A              =[tA EigenvaluesA];

EigenvectorA               = [Eigenvaluestypes.eigenvectors.lambda1 Eigenvaluestypes.eigenvectors.lambda2 Eigenvaluestypes.eigenvectors.lambda3];
EigenvectorA(~NullsEigen.eigA,:)= NaN;
Nulls.Eigenvector.A              =[tA EigenvectorA];

EigenvaluesAs               = Eigenvaluestypes.lambda;
EigenvaluesAs(~NullsEigen.eigAs,:)= NaN;
tAs                         =ttt;
tAs(~NullsEigen.eigAs,1)          =NaN;
Nulls.Eigenvalues.As             =[tAs EigenvaluesAs];

EigenvectorAs               = [Eigenvaluestypes.eigenvectors.lambda1 Eigenvaluestypes.eigenvectors.lambda2 Eigenvaluestypes.eigenvectors.lambda3];
EigenvectorAs(~NullsEigen.eigAs,:)= NaN;
Nulls.Eigenvector.As             =[tAs EigenvectorAs];

EigenvaluesB               = Eigenvaluestypes.lambda;
EigenvaluesB(~NullsEigen.eigB,:)= NaN;
tB                         =ttt;
tB(~NullsEigen.eigB,1)          =NaN;
Nulls.Eigenvalues.B              =[tB EigenvaluesB];

EigenvectorB               = [Eigenvaluestypes.eigenvectors.lambda1 Eigenvaluestypes.eigenvectors.lambda2 Eigenvaluestypes.eigenvectors.lambda3];
EigenvectorB(~NullsEigen.eigB,:)= NaN;
Nulls.Eigenvector.B              =[tB EigenvectorB];

EigenvaluesBs               = Eigenvaluestypes.lambda;
EigenvaluesBs(~NullsEigen.eigBs,:)= NaN;
tBs                        =ttt;
tBs(~NullsEigen.eigBs,1)          =NaN;
Nulls.Eigenvalues.Bs             =[tBs EigenvaluesBs];

EigenvectorBs               = [Eigenvaluestypes.eigenvectors.lambda1 Eigenvaluestypes.eigenvectors.lambda2 Eigenvaluestypes.eigenvectors.lambda3];
EigenvectorBs(~NullsEigen.eigBs,:)= NaN;
Nulls.Eigenvector.Bs              =[tB EigenvectorBs];

Eigenvaluesx               = Eigenvaluestypes.lambda;
Eigenvaluesx(~NullsEigen.eigx,:)= NaN;
tx                         =ttt;
tx(~NullsEigen.eigx,1)          =NaN;
Nulls.Eigenvalues.x              =[tx Eigenvaluesx];

Eigenvectorx               = [Eigenvaluestypes.eigenvectors.lambda1 Eigenvaluestypes.eigenvectors.lambda2 Eigenvaluestypes.eigenvectors.lambda3];
Eigenvectorx(~NullsEigen.eigx,:)= NaN;
Nulls.Eigenvector.x              =[tB Eigenvectorx];

Eigenvalueso               = Eigenvaluestypes.lambda;
Eigenvalueso(~NullsEigen.eigo,:)= NaN;
to                         =ttt;
to(~NullsEigen.eigo,1)          =NaN;
Nulls.Eigenvalues.o              =[to Eigenvalueso];

Eigenvectoro               = [Eigenvaluestypes.eigenvectors.lambda1 Eigenvaluestypes.eigenvectors.lambda2 Eigenvaluestypes.eigenvectors.lambda3];
Eigenvectoro(~NullsEigen.eigB,:)= NaN;
Nulls.Eigenvector.o              =[tB Eigenvectoro];

Eigenvaluesunknown              = Eigenvaluestypes.lambda;
Eigenvaluesunknown(~NullsEigen.unknown,:)= NaN;
tunknown                          =ttt;
tunknown(~NullsEigen.unknown,1)          =NaN;
Nulls.Eigenvalues.unknown        =[tunknown Eigenvaluesunknown];

Eigenvectorunknown               = [Eigenvaluestypes.eigenvectors.lambda1 Eigenvaluestypes.eigenvectors.lambda2 Eigenvaluestypes.eigenvectors.lambda3];
Eigenvectorunknown(~NullsEigen.unknown,:)= NaN;
Nulls.Eigenvector.unknown              =[tB Eigenvectorunknown];

%A
distanceANull                         = dRmin;
distanceANull(~NullsEigen.eigA,2)          = NaN;
distanceANull(~sortdr,2)              = NaN;

positionANull                         = Rn;
positionANull(~NullsEigen.eigA,:)          = NaN;

%B
distanceBNull                         = dRmin;
distanceBNull(~NullsEigen.eigB,2)          = NaN;
distanceBNull(~sortdr,2)              = NaN;

positionBNull                         = Rn;
positionBNull(~NullsEigen.eigB,:)          = NaN;
%X
distanceXNull                         = dRmin;
distanceXNull(~NullsEigen.eigx,2)          = NaN;
distanceXNull(~sortdr,2)              = NaN;

positionXNull                         = Rn;
positionXNull(~NullsEigen.eigx,:)          = NaN;
%Bs
distanceBsNull                        = dRmin;
distanceBsNull(~NullsEigen.eigBs,2)        = NaN;
distanceBsNull(~sortdr,2)             = NaN;

positionBsNull                        = Rn;
positionBsNull(~NullsEigen.eigBs,:)        = NaN;
%As
distanceAsNull                        = dRmin;
distanceAsNull(~NullsEigen.eigAs,2)        = NaN;
distanceAsNull(~sortdr,2)             = NaN;

positionAsNull                        = Rn;
positionAsNull(~NullsEigen.eigAs,:)        = NaN;
%O
distanceONull                         = dRmin;
distanceONull(~NullsEigen.eigo,2)          = NaN;
distanceONull(~sortdr,2)              = NaN;

positionONull                         = Rn;
positionONull(~NullsEigen.eigo,:)          = NaN;
%Unknown type
distanceUnknownNull                   = dRmin;
distanceUnknownNull(~NullsEigen.unknown,2) = NaN;
distanceUnknownNull(~sortdr,2)        = NaN;

positionUnknownNull                         = Rn;
positionUnknownNull(~NullsEigen.unknown,:)       = NaN;

Nulls.NullType.distance.A       = distanceANull;
Nulls.NullType.distance.B       = distanceBNull;
Nulls.NullType.distance.As      = distanceAsNull;
Nulls.NullType.distance.Bs      = distanceBsNull;
Nulls.NullType.distance.unknown = distanceUnknownNull;
Nulls.NullType.distance.x       = distanceXNull;
Nulls.NullType.distance.o       = distanceONull;

Nulls.NullType.position.A       = positionANull;
Nulls.NullType.position.B       = positionBNull;
Nulls.NullType.position.As      = positionAsNull;
Nulls.NullType.position.Bs      = positionBsNull;
Nulls.NullType.position.unknown = positionUnknownNull;
Nulls.NullType.position.x       = positionXNull;
Nulls.NullType.position.o       = positionONull;

Nulls.Requirement.BfieldandEigenvalueslessthanchosenpercentage=constraint;
Nulls.Requirement.DistancewithinSCconfiguration=sortdr;

Nulls.Rnull.Rn1=Rn1;
Nulls.Rnull.Rn2=Rn2;
Nulls.Rnull.Rn3=Rn3;
Nulls.Rnull.Rn4=Rn4;
end

function [NullsEigen,Eigenvalues,Current,constraint,errors]=c_4_null_type(R1,R2,R3,R4,B1,B2,B3,B4,threshold)
%C_4_NULL_TYPE - Determines the type of null
%
%This function determines what type of null is in each data point which
%fulfills the threshold restriction using the nomenclature of Lau and Finn
%(1990) TAJ
%
%   [constraint,Nulls]=C_4_NULL_TYPE(R1,R2,R3,R4,B1,B2,B3,B4,threshold)
%   OUTPUT
%   Nulls is a structure containing logical vectors for each null type
%   constraint is a logical vector that shows true when the threshold
%   restrictions are fulfilled this will be later used to NaN those data
%   points
%   INPUT
%   threshold = tolerance value for accepting data points in percentage.
%   Default=40, threshold=100 means no restriction.
%   B? = the B-field measured at satellite ?: column 1 - time
%   column 2-4 B-field in x,y,z direction
%
% Important: Time tags should be included with the vectors and threshold
% should be given in percentage not deciamalform
%
%See Also C_4_NULL_POSITION
%
%----written by Huishan Fu at BUAA (2014-05-27)----
%----Modified by E.Eriksson----

%Calculates the gradB used in the taylor expansion
gradB=c_4_grad(R1,R2,R3,R4,B1,B2,B3,B4);

%Error in percentage - estimate if linear interpolation is valid to use
[divB,~]=c_4_grad('R?','B?','div');
%jmag=magn_current(:,[1 5]); %fifth column contains the abs value of j (sqrt(j(:,2).^2+j(:,3).^2+j(:,4).^2)) for each time tag
err_4C=irf_multiply(1,real(divB),1,[divB(:,1) real(max(gradB(:,2:end),[],2))],-1); %Essentially does divB/jmag
err_4C(:,2)=abs(err_4C(:,2))*100;
err_4C(:,1)=B1(:,1);

% Only interested in the time intervalls when the two errors (eigenerr and
% curlometer error) is both lower or equal to 40%
eiglamnda=zeros(length(gradB(:,1)),3);
eigenvector1=zeros(length(gradB(:,1)),3);
eigenvector2=zeros(length(gradB(:,1)),3);
eigenvector3=zeros(length(gradB(:,1)),3);
eigVal_err=zeros(length(gradB(:,1)),2);
constraint=false(length(gradB(:,1)),1);
jParallel=zeros(length(gradB(:,1)),2); 
jPerpendicular=zeros(length(gradB(:,1)),2);
%Calculates the parallel and perpendicular current for each null
   [j,~,~ ,~,~,~] = c_4_j(R1,R2,R3,R4,B1,B2,B3,B4);
   temp=irf_abs(j);
for i=1:length(gradB(:,1))
    if all(isnan(gradB(i,2:end)))
        continue
    end
    deltaB_null      = reshape(gradB(i,2:end),3,3);
    [V,D]            = eig(deltaB_null);
    eigenvector1(i,:)=[V(1,1) V(2,1) V(3,1)];
    eigenvector2(i,:)=[V(1,2) V(2,2) V(3,2)];
    eigenvector3(i,:)=[V(1,3) V(2,3) V(3,3)];
    eiglamnda(i,:)   = [D(1,1) D(2,2) D(3,3)];
    eigVal_err(i,2)  = abs(real(D(1,1)+D(2,2)+D(3,3)))/max(abs([real(D(1,1)), real(D(2,2)), real(D(3,3))])) * 100;
    %Calculates the perpendicular and parallel current of the nulls found
    temptime=j(i,1);
   current=temp(i,5);
   if sign(real(D(1,1))) ~= sign(real(D(2,2))) && sign(real(D(1,1))) ~= sign(real(D(3,3)))
       tempEigenVector=irf_abs(eigenvector1(i,:));
       eigenvectors.lambda1=tempEigenVector(1,4);
       jParallel(i,:)=[temptime dot(current,eigenvectors.lambda1)];
       jPerpendicular(i,:)=[temptime current-dot(jParallel(i,2),eigenvectors.lambda1)];
   end
   if sign(real(D(2,2))) ~= sign(real(D(1,1))) && sign(real(D(2,2))) ~= sign(real(D(3,3)))
       tempEigenVector=irf_abs(eigenvector2(i,:));
       eigenvectors.lambda2=tempEigenVector(1,4);
       jParallel(i,:)=[temptime dot(current,eigenvectors.lambda2)];
       jPerpendicular(i,:)=[temptime current-dot(jParallel(i,2),eigenvectors.lambda2)];
   end
   if sign(real(D(3,3))) ~= sign(real(D(1,1))) && sign(real(D(3,3))) ~= sign(real(D(2,2)))
       tempEigenVector=irf_abs(eigenvector3(i,:));
       eigenvectors.lambda3=tempEigenVector(1,4);
       jParallel(i,:)=[temptime dot(current,eigenvectors.lambda3)];
       jPerpendicular(i,:)=[temptime current-dot(jParallel(i,2),eigenvectors.lambda3)];
   end
    if threshold == 100
        constraint(i,1)=true;
    else
        if err_4C(i,2) <= threshold && eigVal_err(i,2) <= threshold
            constraint(i,1)=true;
        else
            constraint(i,1)=false;
        end
    end
end
eiglamnda(~constraint,:)= NaN;
eigenvector1(~constraint,:)= NaN;
eigenvector2(~constraint,:)= NaN;
eigenvector3(~constraint,:)= NaN;
jParallel(~constraint,:)= NaN;
jPerpendicular(~constraint,:)= NaN;
eigVal_err(:,1)=B1(:,1);
eigVal_err(~constraint,:)= NaN;
err_4C(~constraint,:)= NaN;
%Determine the null type by using eigenvalues
%Ideal case
isAllEigenvaluesReal     = abs(max(imag(eiglamnda),[],2)) == 0;
signOfRealpart           = sign(real(eiglamnda));
signOfImaginarypart      = sign(imag(eiglamnda));
nImaginaryNegativeEigenvalues = sum(signOfImaginarypart==-1,2);
nImaginaryPositiveEigenvalues = sum(signOfImaginarypart==+1,2);
nRealNegativeEigenvalues = sum(signOfRealpart==-1,2);
unknowntype              = (nImaginaryNegativeEigenvalues == 3 |sum(signOfImaginarypart,2) ~=0 | nImaginaryPositiveEigenvalues == 3);
aType                    = nRealNegativeEigenvalues == 2 ;
bType                    = nRealNegativeEigenvalues == 1 ;
unknownrealtype          = (nRealNegativeEigenvalues == 3) | (nRealNegativeEigenvalues == 0);
minAbsEigenvalue         = min(abs(eiglamnda),[],2);
twoDType                 = (minAbsEigenvalue == 0);

eigA                     = ( isAllEigenvaluesReal & aType & ~unknowntype & ~unknownrealtype);
eigAs                    = (~isAllEigenvaluesReal & aType & ~unknowntype & ~unknownrealtype);
eigB                     = ( isAllEigenvaluesReal & bType & ~unknowntype & ~unknownrealtype);
eigBs                    = (~isAllEigenvaluesReal & bType & ~unknowntype & ~unknownrealtype);
eigx                     = ( isAllEigenvaluesReal & twoDType & ~unknowntype);
eigo                     = (~isAllEigenvaluesReal & twoDType & ~unknowntype);
unknown                  = (~eigA & ~eigB & ~eigAs & ~eigBs & ~eigx & ~eigo & constraint | unknowntype | unknownrealtype);
NullsEigen.eigA=eigA;
NullsEigen.eigB=eigB;
NullsEigen.eigAs=eigAs;
NullsEigen.eigBs=eigBs;
NullsEigen.eigx=eigx;
NullsEigen.eigo=eigo;
NullsEigen.unknown=unknown;
Eigenvalues.lambda=eiglamnda;
Eigenvalues.eigenvectors.lambda1=eigenvector1;
Eigenvalues.eigenvectors.lambda2=eigenvector2;
Eigenvalues.eigenvectors.lambda3=eigenvector3;
errors.eigval=eigVal_err;
errors.curlometer=err_4C;
Current.jPerpendicular=jPerpendicular;
Current.jParallel=jParallel;
end
