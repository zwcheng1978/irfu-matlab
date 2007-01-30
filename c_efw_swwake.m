function [data, n_corrected] = c_efw_swwake(e,pair,phase_2,plotflag)
%C_EFW_SWWAKE  Correct raw EFW E data for wake in the solar wind
%
% [data, n_spins_corrected] = c_efw_swwake(e,pair,phase_2 [,plotflag])
%
% Wakes are identified by max derivative
%
% plotflag=0, no plotting of debug 
%
% $Id$

% Copyright 2007 Yuri Khotyaintsev

error(nargchk(3,4,nargin))
if nargin==4
	if isnumeric(plotflag)
		if plotflag~=0, plotflag = 1; end
	else plotflag = 1;
	end
else plotflag = 0;
end

n_corrected = 0;

if pair~=12 && pair~=32 && pair~=34, error('PAIR must be one of: 12, 32, 34'), end
if size(phase_2,1)<2, error('not enough points in phase_2'), end

if pair==32, error('PAIR 32 is not implemented yet'), end

% N_EMPTY .75 means that we use only spins with more then 75% points.
N_EMPTY = .9; 
MAX_SPIN_PERIOD = 4.3;
WAKE_MAX_HALFWIDTH = 45; % degrees
WAKE_MIN_AMPLITUDE = 0.1; % mV/m

data = e;

i0 = find(phase_2(:,2)==0);
if isempty(i0), irf_log('proc','empty phase'), return, end
if length(i0)<5, irf_log('proc','not enough spins for correction'), return, end

% Hack to add phase=360 at the end and phase=0 at the start of the interval
% in order to try to correct these incomplete spins
if phase_2(end,2)~=360 && ...
		phase_2(i0(end),1)-phase_2(i0(end-1),1)<MAX_SPIN_PERIOD
	phase_2 = [phase_2; phase_2(i0(end),1)*2-phase_2(i0(end-1),1) 0];
end
if phase_2(1,2)~=0 && phase_2(i0(2),1)-phase_2(i0(1),1)<MAX_SPIN_PERIOD
	phase_2 = [phase_2(i0(1),1)*2-phase_2(i0(2),1) 0; phase_2];
	i0 = find(phase_2(:,2)==0);
end

n_spins = length(i0);
NPOINTS = 361;
tt = zeros(NPOINTS,length(i0));
ttime = tt;
sf = c_efw_fsample(e,'hx');
iok = [];

for in = 1:n_spins
	ts = phase_2(i0(in),1);
	i360 = find( phase_2(:,1)>ts & phase_2(:,1)<ts+MAX_SPIN_PERIOD & ...
		phase_2(:,2)==360);
	if isempty(i360)
		irf_log('proc',['gap in phase at ' epoch2iso(ts,1)])
		te = ts + 4.0;
		empty = 1;
	else
		if length(i360)~=1, error('bogus phase'), end
		te = phase_2(i360,1);
		empty = 0;
	end
	
	ttime(:,in) = (ts + (0:1:360) *(te-ts)/360.0)';
	if empty, tt(:,in) = NaN;
	else
		eind = find((e(:,1) > ts-0.15) & (e(:,1) < te+0.15));
		eind(isnan(e(eind))) = [];
		% Check for data gaps inside one spin.
		if sf>0 && length(eind)<N_EMPTY*MAX_SPIN_PERIOD*sf
			irf_log('proc',['data gap at ' epoch2iso(ts,1)])
			tt(:,in) = NaN;
		else
			if sf==450, dtmp = c_resamp(e(eind,:), ttime(:,in));
			else dtmp = c_resamp(e(eind,:), ttime(:,in), 'spline');
			end
			tt(:,in) = dtmp(:,2);
		end
		
		% Identify spins for which we attempt to correct wake
		if (in>=6 && sum(isnan(tt(1,in - (0:1:5) )))<=1)
			% We allow max one data gap within 6 spins
			iok = [iok in-3];
		end
		if (in==n_spins || n_spins==5) && ~any(isnan(tt(1,in - (1:4) )))
			% Special case when we have only 5 spins or
			% it is the last spin we are working with
			iok = [iok in-2];
		end
	end
end

if isempty(iok)
	irf_log('proc','not enough spins for correction')
	return
end

i1 = [];

for in = iok
	ts = ttime(1,in);
	if in==n_spins-2
		% Last interval
		idx = -2:1:2;
	else
		% Check for a data gap
		nans = isnan(tt(1,in + (-2:1:3) ));
		if any(nans)
			idx = 1:6;
			idx = idx(xor(idx,nans)) -3;
		else idx = -2:1:2;
		end
	end
	
	av12 = mean(tt(:, in + (idx) ),2);

	% Identify wakes by max derivative
	d12 = diff(diff(av12));
	d12 = [d12(1); d12; d12(end)];
	if plotflag
		d12_tmp = d12; % save for plotting
	end
	d12 = w_ave(d12);
	
	% If i1 is coming not from the previous spin
	if isempty(i1) || (in~=iok(1) && in-1~=iok(find(iok==in)-1))
		i1 = 1:1:361;
	end
	
	% We expect for the second maximum to be 180 degrees from the first
	if max(d12(i1))>abs(min(d12(i1)))
		ind1 = find(d12 == max(d12(i1)))+1;
		i1 = mod( ind1-WAKE_MAX_HALFWIDTH+180:ind1+WAKE_MAX_HALFWIDTH+180,...
			NPOINTS) +1;
		ind2 = find(d12 == min(d12( i1 )))+1;
	else
		ind1 = find(d12 == min(d12(i1)))+1;
		i1 = mod( ind1-WAKE_MAX_HALFWIDTH+180:ind1+WAKE_MAX_HALFWIDTH+180,...
			NPOINTS) +1;
		ind2 = find(d12 == max(d12( i1 )))+1;
	end
	
	i1 = mod( ind1-WAKE_MAX_HALFWIDTH:ind1+WAKE_MAX_HALFWIDTH , NPOINTS) +1;
	i2 = mod( ind2-WAKE_MAX_HALFWIDTH:ind2+WAKE_MAX_HALFWIDTH , NPOINTS) +1;
	
	dav = (d12(i1)-d12(i2))/2;
	cdav = cumsum(dav);
	cdav = cdav - mean(cdav);
	ccdav = cumsum(cdav);
	if max(abs(ccdav))< WAKE_MIN_AMPLITUDE
		irf_log('proc',['wake is too small at ' epoch2iso(ts,1)])
		i1 = [];
		continue
	end
	wake = zeros(size(av12));
	wake(mod(i1,NPOINTS) +1) = ccdav;
	wake(mod(i2,NPOINTS) +1) = -ccdav;
	
	if plotflag
		subplot(4,1,1)
		plot(ttime(:,in)-ts, tt(:, in + (-2:1:2) ), 'g',...
			ttime(:,in)-ts, av12, 'k',...
			ttime(ind1+1,in)*[1 1]-ts, [-2 2], 'r',...
			ttime(ind2+1,in)*[1 1]-ts, [-2 2], 'r',...
			ttime(:,in)-ts, av12-wake,'r');
		ylabel('E12 [mV/m]');
		add_timeaxis(gca,ts); xlabel('');

		subplot(4,1,2)
		plot(ttime(:,in)-ts,d12_tmp,'g',ttime(:,in)-ts, d12,'b');
		ylabel(['D2(E' num2str(pair) ') [mV/m]']);
		add_timeaxis(gca,ts); xlabel('');

		subplot(4,1,3)
		plot(ttime(:,in)-ts, wake)
		ylabel('Wake [mV/m]');
		add_timeaxis(gca,ts);
	end
	
	ind = find(e(:,1)>=ttime(1,in) & e(:,1)<ttime(end,in));
	wake_e = c_resamp([ttime(:,in) wake], e(ind,1));
	data(ind,2) = data(ind,2) - wake_e(:,2);
	n_corrected = n_corrected + 1;
	
	if plotflag
		subplot(4,1,4)
		ts = ttime(1,in);
		te = ttime(end,in);
		irf_plot({e(ind,:),data(ind,:)},'comp')
		ylabel(['E' num2str(pair) ' [mV/m]']);
	end
	
	cox = [];
	if in==iok(1) || (in~=iok(1) && in-1~=iok(find(iok==in)-1))
		% If the previous spin was not corrected
		% we correct it here
		if in==iok(1) && in>3
			% We try to correct the spin at the start of the entire 
			% intrval which usually has a data gap 
			cox = - (1:3);
			n_corrected = n_corrected + 3;
		else
			cox = - (1:2);
			n_corrected = n_corrected + 2;
		end
	end
	if in==iok(end) || (in~=iok(end) && in+1~=iok(find(iok==in)+1))
		if ~isempty(cox)
			irf_log('proc',['single interval at ' epoch2iso(ts,1)])
		else
			irf_log('proc',['stop   interval at ' epoch2iso(ts,1)])
		end
		if in==iok(end) && n_spins-in>=3
			% We try to correct the spin at the end of the entire 
			% intrval which usually has a data gap 
			cox = [cox 1:3];
			n_corrected = n_corrected + 3;
		else
			cox = [cox 1:2];
			n_corrected = n_corrected + 2;
		end
	elseif ~isempty(cox)
		irf_log('proc',['start  interval at ' epoch2iso(ts,1)])
	end
	if ~isempty(cox)
		for cx = cox
			ind = find(e(:,1)>=ttime(1,in+cx) & e(:,1)<ttime(end,in+cx));
			wake_e = c_resamp([ttime(:,in+cx) wake], e(ind,1));
			data(ind,2) = data(ind,2) - wake_e(:,2);
			%irf_log('proc',['correcting spin: ' num2str(cx)])
			if plotflag
				hold on
				irf_plot({e(ind,:),data(ind,:)},'comp')
				ylabel(['E' num2str(pair) ' [mV/m]']);
				hold off
				if ttime(1,in+cx)<ts, ts = ttime(1,in+cx); end
				if ttime(end,in+cx)>te, te = ttime(end,in+cx); end
				irf_zoom([ts te],'x',gca)
			end
		end
	end
	
end

irf_log('proc',['corrected ' num2str(n_corrected) ' out of ' ...
	num2str(n_spins) ' spins'])

function av = w_ave(x)
% Weinghted average
NPOINTS = 361;
m = [.07 .15 .18 .2 .18 .15 .07]';
av = zeros(size(x));
for j=1:length(x)
	ii = j + (-3:1:3);
	if j<=3, ii(ii<1) = ii(ii<1) +NPOINTS; end
	if j>length(x)-3, ii(ii>NPOINTS) = ii(ii>NPOINTS) -NPOINTS; end
	av(j) = sum(x(ii).*m);
end