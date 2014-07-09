function [constraint,Nulls]=c_4_null_type(R1,R2,R3,R4,B1,B2,B3,B4,threshold)
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
%   restrictions wasn't fulfilled this will be later used to NaN those data
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
[j,divB] = c_4_j('R?','B?');
magn_current=irf_abs(j);
jmag=magn_current(:,[1 5]); %fifth column contains the abs value of j (sqrt(j(:,2).^2+j(:,3).^2+j(:,4).^2)) for each time tag
err_4C=irf_multiply(1,divB,1,jmag,-1); %Essentially does divB/jmag
err_4C(:,2)=abs(err_4C(:,2))*100;

% Only interested in the time intervalls when the two errors (eigenerr and
% curlometer error) is both lower or equal to 40%
eigvec=zeros(length(gradB(:,1)),3);
eigVal_err=zeros(length(gradB(:,1)),2);
constraint=false(length(gradB(:,1)),1);
for i=1:length(gradB(:,1))
    if all(isnan(gradB(i,2:end)))
        continue
    end
    deltaB_null     = reshape(gradB(i,2:end),3,3);
    D               = eig(deltaB_null);
    eigvec(i,:)     = [D(1,1) D(2,1) D(3,1)];
    eigVal_err(i,2) = abs(real(D(1,1)+D(2,1)+D(3,1)))/max(abs([real(D(1,1)), real(D(2,1)), real(D(3,1))])) * 100;
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
eigvec(~constraint,:)= NaN;

%Determine the null type by using eigenvalues
%Ideal case
isAllEigenvaluesReal     = max(imag(eigvec),[],2) == 0;
signOfRealpart           = sign(real(eigvec));
nRealNegativeEigenvalues = sum(signOfRealpart==-1,2);
aType                    = (nRealNegativeEigenvalues == 2);
bType                    = (nRealNegativeEigenvalues == 1);
minAbsEigenvalue         = min(abs(eigvec),[],2);
twoDType                 = (minAbsEigenvalue == 0);

eigA                     = ( isAllEigenvaluesReal & aType);
eigAs                    = (~isAllEigenvaluesReal & aType);
eigB                     = ( isAllEigenvaluesReal & bType);
eigBs                    = (~isAllEigenvaluesReal & bType);
eigx                     = ( isAllEigenvaluesReal & twoDType);
eigo                     = (~isAllEigenvaluesReal & twoDType);
unknown                  = (~eigA & ~eigB & ~eigAs & ~eigBs & ~eigx & ~eigo & constraint);

Nulls.eigA=eigA;
Nulls.eigB=eigB;
Nulls.eigAs=eigAs;
Nulls.eigBs=eigBs;
Nulls.eigx=eigx;
Nulls.eigo=eigo;
Nulls.unknown=unknown;
end


