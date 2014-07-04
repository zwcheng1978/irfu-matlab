function index=c_4_poincare_index(B1,B2,B3,B4)
%%C_4_POINCARE_INDEX - Calculates the poincare index using B-field measurements
%from four S/C's in a tetrahedron formation. 
%
%This function calculates the poincar� index from B-fields measurements. If index is a nonzero number
%than there is at least one null point within the volume.
%
%   index=poincare_index_3D(B1,B2,B3,B4)
%   index = [time index]
%   B? is the B-field measured at satellite ?
%
% Important: This method only works in the 3D case. That is when none of
% the eigenvalues in deltaB vanish. The time tags needs to be included in each vector or the program won't
% work
% See Also IRF.SOLIDANGLE
% Reference: Greene 1992 JCP (98) p.194-198

%%------written by Yunhui Hu, Dec.20.2006 in Wuhan------------
%------modified by Shiyong Huang, May.17.2012 in IRF---------
%------modified by E.Eriksson in IRFU------------

n=size(B1,2); 
if n<4
error('Time tag must be included in each input vector. Please do so and try again.')
end
if nargin==0
    help c_4_null_position;
    return;
elseif nargin < 4
    error('Too few input values. See usage: help poincare_index_3D')
elseif nargin>4
    error('Too many input values. See usage: help poincare_index_3D')
end

%Each vector contains time tags so all vectors needs to be resampled to
%establish synchronisation between the S/C's
if  isempty(isnan(B1(:,1)))
B2=irf_resamp(B2,B1);
B3=irf_resamp(B3,B1);
B4=irf_resamp(B4,B1);
end

%Renaming of vectors for simplicty sake and removing time for calculations
%using B?
time=B1(:,1); %time= first column of SC1
B1=B1(:,2:4); %vec1= 2 to 4th column (Bx,By,Bz) for s/c 1 (if that's the one you placed in position 1)
B2=B2(:,2:4); %vec2= 2 to 4th column (Bx,By,Bz)
B3=B3(:,2:4); %vec3= 2 to 4th column (Bx,By,Bz)
B4=B4(:,2:4); %vec4= 2 to 4th column (Bx,By,Bz)

lx=size(B1,1); %lx becomes the number of rows in B1

%Create the zero matrix that will be used to map each point from xyz to
%By,Bx,Bz space

Map_sc1=zeros(lx,3); 
Map_sc2=zeros(lx,3); 
Map_sc3=zeros(lx,3);
Map_sc4=zeros(lx,3);
% map the points from x,y,z to magnetic three-dimensional field space (Bx,By,Bz instead of
% x,y,z) (ref. Greene, J.M. 1990)
% A null point in configuration space (xyz satellites) corresponds to the
% origin in M space.

%Mapping of the B fields for each s/c but as a unit vector (length=1)
%Calculate the length of each vector
   vl1= sqrt(dot(B1,B1,2)); %norm (length of vector)=sqrt(dot(vec1,vec1,2)).
   %dot(vec1,vec1,2) treats each row as a vector in the matrix so A=dot(vec1,vec1,2)
   % would give A(1,:) dot product of vec1(1,:) and vec1(1,:)
   vl2= sqrt(dot(B2,B2,2)); 
   vl3= sqrt(dot(B3,B3,2)); 
   vl4= sqrt(dot(B4,B4,2)); 
   
   %%Unit vectors is used in solid angle so we need divide the vector with
   %%its norm (magnitude of the vector/length)
   for i=1:3    %Each column on Map_sc1 is given unit vector (first column in 3D divided by the length of vector in the 
       % direction of each s/c)
    Map_sc1(:,i)=B1(:,i)./vl1; %Needs to use the for loop so that matrix dimensions agrees
    Map_sc2(:,i)=B2(:,i)./vl2; 
    Map_sc3(:,i)=B3(:,i)./vl3; 
    Map_sc4(:,i)=B4(:,i)./vl4; 
   end
    %Calculate the solid angle of an unit sphere that has taken the sign into account to give the number of nulls for each triangles with the direction going
    %counterclockwise seen from origin which determines which sc you choose
    %as point A,B,C see solidangle for more details
    area1=irf.solidangle(Map_sc1,Map_sc2,Map_sc3); %In a tetrahedron you have four triangles so you need the solid angle for each
    area2=irf.solidangle(Map_sc1,Map_sc4,Map_sc2);
    area3=irf.solidangle(Map_sc1,Map_sc3,Map_sc4);
    area4=irf.solidangle(Map_sc2,Map_sc4,Map_sc3);
    %The Poincar? index is the total area (including each sign) divided by
    %4pi
    index=(area1+area2+area3+area4)/(4*pi); 
    
    %add time tags
    index=[time index]; %If the indices is a nonzero then there is at least one null point within the volume of interest
end