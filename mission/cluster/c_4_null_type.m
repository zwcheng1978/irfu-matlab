function [constraint,Nulls]=c_4_null_type(threshold,B1,B2,B3,B4,R1,R2,R3,R4)
%C_4_NULL_TYPE - Determines the type of null
%
%This function determines what type of null is in each data point which
%fulfills the threshold restriction using the nomenclature of Lau and Finn
%(1990) TAJ
%
%   [constraint,Nulls]=c_4_null_type(threshold,B1,B2,B3,B4,R1,R2,R3,R4)
%   Nulls is a structure containing logical vectors for each null type
%   constraint is a logical vector that shows true when the threshold
%   restrictions wasn't fulfilled this will be later used to NaN those data
%   points
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
constraint=true(length(gradB(:,1)),1);
for i=1:length(gradB(:,1))
    deltaB_null=reshape(gradB(i,2:end),3,3);
    D = eig(deltaB_null);
    eigvec(i,:)=[D(1,1) D(2,1) D(3,1)];
    eigVal_err(i,2)=abs(real(D(1,1)+D(2,1)+D(3,1)))/max(abs([real(D(1,1)), real(D(2,1)), real(D(3,1))])) * 100;
    if threshold==100
        constraint(i,1)=false;
    else
    if err_4C(i,2) <= threshold
        if eigVal_err(i,2) <= threshold
            constraint(i,1)=false;
        else
            constraint(i,1)=true;
            
        end
    else
        constraint(i,1)=true;
    end
    end
end
t=gradB(:,1);
t(constraint,:)=NaN;
eigvec(constraint,:)= NaN;

%Determine the null type by using eigenvalues
eigA=true(length(t(:,1)),1);
eigAs=true(length(t(:,1)),1);
eigB=true(length(t(:,1)),1);
eigBs=true(length(t(:,1)),1);
eigx=true(length(t(:,1)),1);
eigo=true(length(t(:,1)),1);
unknown=true(length(t(:,1)),1);

%Ideal case
for i=1:length(t(:,1))
    if max(abs([imag(eigvec(i,1)) imag(eigvec(i,2)) imag(eigvec(i,3))])) == 0  %Checks that all values are real (Type A or B)
        if length(find([eigvec(i,1) eigvec(i,2) eigvec(i,3)]<0)) == 2  %If two of the real values are negative Type A
            eigA(i,1)=false;
        else
            if length(find([eigvec(i,1) eigvec(i,2) eigvec(i,3)]<0)) == 1 %Checks if only one real values is positive Type B
                eigB(i,1)=false;
            else
                unknown(i,1)=false;
            end
        end
        if min(abs([eigvec(i,1) eigvec(i,2) eigvec(i,3)]))==0  %Checks if one of the eigenvalues are zero therefore a 2D null
            eigx(i,1)=false;
        end
        %If not all of the eigenvalues are real then we have either a Bs,As
        %or O null type
    else
        if length(find([real(eigvec(i,1)) real(eigvec(i,2)) real(eigvec(i,3))]>0)) == 2    %Checks how many positive real eigenvalues there is a Bs contains two positive and one neg.
            %As contains two neg and one pos and O only contains imaginary
            %values
            eigBs(i,1)=false;
        else
            if length(find([real(eigvec(i,1)) real(eigvec(i,2)) real(eigvec(i,3))]>0)) == 1    %Are there only 1 positive real eigenvalue then the type is As
                eigAs(i,1)=false;
            else
                unknown(i,1)=false;
            end
        end
        if max(abs([real(eigvec(i,1)) real(eigvec(i,2)) real(eigvec(i,3))]))==0    %If there is no real eigenvalue then you have a 0 eigenvalue
            eigo(i,1)=false;
        end
    end
    Nulls.eigA=eigA;
    Nulls.eigB=eigB;
    Nulls.eigAs=eigAs;
    Nulls.eigBs=eigBs;
    Nulls.eigx=eigx;
    Nulls.eigo=eigo;
    Nulls.unknown=unknown;
end

   
