function ok=test(varargin)
% IRF.TEST test different irfu-matlab routines
%
%   IRF.TEST show list of possible tests
%   IRF.TEST('full')  - make full test
%	IRF.TEST(testname) - make test 'testname'
%
%	Example:
%   IRF.TEST('irf_time')  - test irf_time
%

% $Id$

testList={'irf_time','epoch','c_4','coord_sys'};

if nargin == 0,
	disp('Possible tests');
	disp('--------------------------');
	for j=1:numel(testList),
		disp([num2str(j) '. ' testList{j}]);
	end
	disp('--------------------------');
	disp('Execute: irf.test(''testname'') or irf.test(number) or irf.test(''all'')');
	return
else
	testName = varargin{1};
	if ischar(testName),
		isTest = strcmpi(testList,testName);
		if any(isTest)
			testNumber = find(isTest);
		elseif strcmpi('all',testName)
			testNumber=1:numel(testList);
		else
			disp(['!!! Test ''' testName ''' not in list!!!']);
			irf.test;
			return;
		end
	elseif isnumeric(testName) %assume given test number
		testNumber=testName;
	else
		irf_log('fcal','Unknown test');
		return;		
	end
end
nTests = numel(testNumber);
okArray = false(nTests,1);
for j=1:nTests
	iTest=testNumber(j);
	disp(' ');
	disp('*********************');
	if numel(testNumber)>1,
		disp(['TEST ' num2str(j) '/' num2str(nTests) ': ' testList{iTest}]);
	else
		disp(['TEST: ' testList{iTest}]);
	end
	disp('*********************');
	tic;
	functionName = ['test_' testList{iTest}];
	okArray(j)=eval(functionName);
	tElapsed = toc;
	disp(['Test ended. Elapsed time: ' num2str(tElapsed) 's.']);
end

if all(okArray),
	disp('ALL TESTS PASSED SUCCESFULL!');
else
	for j=find(~okArray)
		disp(' ');
		disp(['Test FAILED! ' testList{testNumber(j)}]);
	end
end

if nargout == 1, 
	ok = okArray;
end
end

% Tests
function okTest = test_irf_time
try
	okTest = 1;
	% generate vector with 10000 times during last 500 years
	a=rand(1000,1);
	tDateArray = now - 365*500*a;
	s1=irf_time(tDateArray,'date2iso');
	t=irf_time(s1,'iso2epoch');
	s2=irf_time(t,'epoch2iso');
	ok = strcmp(s1,s2);
	plusminus(ok); disp('1000 random iso>epoch>iso ');	
	okTest = okTest * ok;
	% next test
	tint=[t t+a*1000];
	s1=irf_time(tint,'tint2iso');
	tt=irf_time(s1,'iso2tint');
	s2=irf_time(tt,'tint2iso');
	ok = strcmp(s1,s2);
	plusminus(ok); disp('1000 random iso>tint>iso ');	
	okTest = okTest * ok;
	% next test. different iso formats should be recognized
	s1=irf_time(tDateArray,'date2iso');
	t=irf_time(s1,'iso2epoch');
	t1=irf_time(s1(:,1:end-1),'iso2epoch');
	t2=irf_time(reshape(strrep(s1(:)','T',' '),size(s1)),'iso2epoch');
	t3=irf_time(reshape(strrep(s1(:)','Z',' '),size(s1)),'iso2epoch');
	if all(t==t1) && all(t==t2) && all(t==t3)
		ok=1;
	else
		ok=0;
	end
	plusminus(ok); disp('1000 random iso (4 formats) > tint');	
	okTest = okTest * ok;
	% next test
catch
	okTest=false;
end
end
function okTest = test_c_4
try
	okTest = true; % default for all test
	ok     = true; % default for subtest
	
	%% SUBTEST: position of tetrahedron, base plane in xy
	R1n=[1 0 0];
	R2n=[cos(1/3*2*pi) sin(1/3*2*pi) 0];
	R3n=[cos(2/3*2*pi) sin(2/3*2*pi) 0];
	R4n=[0 0 sqrt(norm(R2n-R1n)^2-1)];
	zMassCenter=R4n(3)/4;
	c_eval('R?n(3)=R?n(3)-zMassCenter;');
	c_eval('R?=R?n;');
	
	% construction of field
	Bconst	=@(x) [0 0 1]; %
	
	% Cylindrical field with symmetry axis along z and center at rc
	% such field gives curvature equal to 1/R where R is distance to
	% the cylinder symmetry axis
	Bcircle	=@(x,rc) [(x(:,2)-rc(2))./sqrt((x(:,1)-rc(1)).^2 +(x(:,2)-rc(2)).^2) ...
		-(x(:,1)-rc(1))./sqrt((x(:,1)-rc(1)).^2 +(x(:,2)-rc(2)).^2) 0];
	rc=[10 10 0]; % cylinder symmetry axis location, B field is clockwise
	%disp('Test curvature')
	%disp([' should be: [' num2str(1/norm(rc)^2*rc,'%5.2f') ']']);
	c_eval('B?=Bcircle(R?,rc);');
	curv=c_4_grad('R?','B?','curvature');
	%disp(['        is: [' num2str(curv,'%5.2f') ']'])
	if norm(1/norm(rc)^2*rc - curv) > 1e-1
		disp('!!!! Curvature failed!!!!')
		disp([' should be: [' num2str(1/norm(rc)^2*rc,'%5.2f') ']']);
		disp(['        is: [' num2str(curv,'%5.2f') ']'])
		disp(['     error: ' num2str(norm(1/norm(rc)^2*rc - curv))])
		ok = false;
	end
	plusminus(ok); disp('curvature');	
	okTest = okTest * ok;
	
	%% SUBTEST: Constant current in Z direction
	jz=10;
	Bjz = @(x) [0 jz*x(1) 0];
	c_eval('B?=Bjz(R?);');
	testCurl=c_4_grad('R?','B?','curl');
	if norm([0 0 jz] - testCurl) > 1e-10
		disp('!!!! Curl failed!!!!')
		disp([' should be: [0 0 ' num2str(jz) ']']);
		disp(['        is: [' num2str(testCurl,'%6.2f') ']'])
		disp(['     error: ' num2str(norm([0 0 jz] - testCurl))])
		ok = false;
	end
	plusminus(ok); disp('curl');	
	okTest = okTest * ok;
	
	% Test drift gradient
	% Fitzpatrick book page 30
	% V=mv^2/2/e B x gradB / B^3

catch
	okTest=false;
end
end
function okTest = test_epoch
try 
	okTest		= true; % default for all test
	ntests=1e5;
	disp('Testing ISDAT epoch conversion tools.')
	disp([' Using ' num2str(ntests) ' times for random tests'])
	
	%% SUBTEST: roundoff error
	okSubtest	= true; % default for subtest
	tol=iso2epoch('2020-01-01T00:00:00Z')*eps;
	t=1272094620.1;  % Known to cause problems on earlier Matlab versions
	d=fromepoch(t);
	if ~isequal(d(1:5),[2010 04 24 07 37]) || abs(d(6)-0.1)>tol
		okSubtest = false;
	end
	plusminus(okSubtest); disp('known earlier roundoff error');	
	okTest = okTest * okSubtest;
	%% SUBTEST: roundoff errors from random times
	okSubtest	= true; % default for subtest
	t=rand(1,ntests)*iso2epoch('2020-01-01T00:00:00Z');
	d=fromepoch(t);
	if any(d(:,6)>60)
		okSubtest = false;
	end
	plusminus(okSubtest); disp('random roundoff errors');	
	okTest = okTest * okSubtest;
	
	%% SUBTEST: roundoff errors from times near second boundaries
	okSubtest	= true; % default for subtest
	t=rand(1,ntests)*iso2epoch('2020-01-01T00:00:00Z');
	t=fix(t)+rand(1,ntests)*0.0001;
	d=fromepoch(t);
	if any(d(:,6)>60)
		okSubtest = false;
	end
	plusminus(okSubtest); disp('roundoff errors near second boundaries');	
	okTest = okTest * okSubtest;
	%% SUBTEST: epoch2iso/iso2epoch round trip
	okSubtest	= true; % default for subtest
	year=fix(rand(ntests,1)*(2020-1970)+1970);
	month=fix(rand(ntests,1)*12)+1;
	day=fix(rand(ntests,1).*(eomday(year,month)))+1;
	hour=fix(rand(ntests,1)*24);
	minute=fix(rand(ntests,1)*60);
	second=rand(ntests,1);
	t=[year month day hour minute second];
	ISO =sprintf('%04d-%02d-%02dT%02d:%02d:%09.6fZ',t');
	ISO =reshape(ISO,27,ntests)';
	ISO2=epoch2iso(iso2epoch(ISO));
	if strcmp(ISO,ISO2)~=1
		okSubtest = false;
	end
	plusminus(okSubtest); disp('epoch2iso/iso2epoch round trip');	
	okTest = okTest * okSubtest;
	if ~okSubtest
		cnt=0;
		for i=1:ntests
			if strcmp(ISO(i,:),ISO2(i,:))~=1
				disp(['  Failed for time ' ISO(i,:) ' --> ' ISO2(i,:)])
				cnt=cnt+1;
			end
			if cnt > 5
				disp('  ... and possibly more times (not listing all failures).')
				break
			end
		end
	end
	
	%% SUBTEST: leap seconds
	okSubtest	= true; % default for subtest
	stepdates = [...
		'Jan 6 1980'
		'Jul 1 1981'
		'Jul 1 1982'
		'Jul 1 1983'
		'Jul 1 1985'
		'Jan 1 1988'
		'Jan 1 1990'
		'Jan 1 1991'
		'Jul 1 1992'
		'Jul 1 1993'
		'Jul 1 1994'
		'Jan 1 1996'
		'Jul 1 1997'
		'Jan 1 1999'
		'Jan 1 2006'
		'Jan 1 2009'];
	stepdates = datenum(stepdates)';
	ISOtime=epoch2iso(date2epoch(stepdates)-0.5);
	for i=1:length(stepdates)
		leap_s=0;
		if strcmp(ISOtime(1,18:19),'60')
			okSubtest = false;
			leap_s=1;
			
			break;
		end
	end
	plusminus(okSubtest); disp('leap seconds');	
	okTest = okTest * okSubtest;
	if leap_s==0,disp('  Neither Matlab nor ISDAT uses leap seconds.'),end
	if leap_s==1 && ~okSubtest, disp('ISDAT does not use leap seconds.'); end


catch
	okTest = false;
end
end
function okTest = test_coord_sys
try 
	okTest		= true; % default for all test	
	%% SUBTEST: conversion between geo/gei/gse/gsm/sm/mag
	okSubtest	= true; % default for subtest
	%
	% Here comes the subtest1
	%
	coordsysList = {'gei','geo','gse','gsm','sm','mag'};
	% generate vector with 1000 times during last 50 years
	iTimes = 100;
	a=rand(iTimes,1);a=a(:);
	tDateArray = now - 365*50*a;
	t=irf_time(tDateArray,'datenum2epoch');
	for iT = 1:numel(t)
		iCoord = [randi(numel(coordsysList),1,3) 0];
		iCoord(end)=iCoord(1);
		vecStart = rand(1,3)*rand(1)*1e4;
		vec=[t(iT) vecStart];
		for jCoord = 1:numel(iCoord)-1,
			transformString = [coordsysList{iCoord(jCoord)} '>' coordsysList{iCoord(jCoord+1)}];
			vec = irf.geocentric_coordinate_transformation(vec,transformString);
		end
		if abs(vec(2:4)-vecStart)>1e-10,
			okSubtest = false;
			disp(['failed time: ' irf_time(t(iT),'iso') ]);
			disp(['failed conversion: ' coordsysList(iCoord) ]);
			disp(['failed start vector: ' num2str(vecStart,'%9.2e')]);
			disp(['failed end vector: ' num2str(vec(2:4),'%9.2e')]);
			disp(['eps: ' num2str(abs(vec(2:4)-vecStart))]);
			break;
		end
	end
	plusminus(okSubtest); disp([num2str(iTimes) ' random cyclic transformations gei/geo/gse/gsm/sm/mag']);	
	okTest = okTest * okSubtest;	
catch
	okTest = false;
end
end
function okTest = template_test
try 
	okTest		= true; % default for all test
	okSubtest	= true; % default for subtest
	
	%% SUBTEST: description of subtest1
	okSubtest	= true; % default for subtest
	%
	% Here comes the subtest1
	%
	plusminus(okSubtest); disp('subtest1 text');	
	okTest = okTest * okSubtest;
	%% SUBTEST: description of test2
	okSubtest	= true; % default for subtest
	%
	% Here comes the subtest2
	%
	plusminus(okSubtest); disp('subtest2 text');	
	okTest = okTest * okSubtest;
	
catch
	okTest = false;
end
end

% Subfunctions
function plusminus(ok)
if ok,
	fprintf('+ ');
else
	fprintf('- ');
end
end

