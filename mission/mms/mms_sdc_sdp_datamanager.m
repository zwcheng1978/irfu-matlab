function [ varOut ] = mms_sdc_sdp_datamanager( param, dataObj )
% mms_sdc_sdp_datamanager will store and retrive data for
%  	mms_sdc_sdp_datamanager( dataType, dataObj ) will store
%  	appropriate data variables related to dataType in the global variable
%  	DATAC for later retrival.
%  
%   [varOut] = mms_sdc_sdp_datamanager( variable ) will return the variable
%   requested to varOut, if no such variable has been created already it
%   will try and calculate it using the stored data.
%   
%  	Example:
%
%  		mms_sdc_sdp_datamanager('dce',dceDataObj)
%  		mms_sdc_sdp_datamanager('phase')
%  
%   	See also DATAOBJ, MMS_SDC_SDP_CDF_IN_PROCESS.

narginchk(1,2); % One argument is simply retreive variable, Two arguments
% store "dataObj" as "dataType".

global MMS_CONST, if isempty(MMS_CONST), MMS_CONST = mms_constants(); end
global DATAC; % Here we store read data.

if nargout, varOut = []; end

if ~ischar(param),
  errStr = 'PARAM must be a string';
  irf.log('critical', errStr); error(errStr);
end

if strcmpi(param, 'init')
  % Initialize
  if nargin==1
    errStr = 'INIT requires second input argument';
    irf.log('critical', errStr); error(errStr);
  elseif ~isstruct(dataObj)
    errStr = 'Second input argument for INIT must be a structure';
    irf.log('critical', errStr); error(errStr);
  elseif ~isfield(dataObj,'scId') || ~isnumeric(dataObj.scId) || ...
      isempty(intersect(dataObj.scId, MMS_CONST.MMSids))
    errStr = 'Invalid input for init_struct.scId';
    irf.log('critical', errStr); error(errStr);
  end
  DATAC.scId = dataObj.scId;
  if ~isfield(dataObj,'tmMode')
    DATAC.tmMode = 1;
    irf.log('warining',['init_struct.tmMode not specified, defaulting to '''...
      MMS_CONST.TmModes{DATAC.tmMode} ''''])
  elseif ~isnumeric(dataObj.tmMode) || ...
      isempty(intersect(dataObj.tmMode, 1:numel(MMS_CONST.TmModes)))
    errStr = 'Invalid input for init_struct.tmMode';
    irf.log('critical', errStr); error(errStr);
  else DATAC.tmMode = dataObj.tmMode;
  end
  if ~isfield(dataObj,'procId')
    DATAC.procId = 1;
    irf.log('warining',['init_struct.procId not specified, defaulting to '''...
      MMS_CONST.SDCProcs{DATAC.procId} ''''])
  elseif ~isnumeric(dataObj.procId) || ...
      isempty(intersect(dataObj.procId, 1:numel(MMS_CONST.SDCProcs)))
    errStr = 'Invalid input for init_struct.procId';
    irf.log('critical', errStr); error(errStr);
  else DATAC.procId = dataObj.procId;
  end
  DATAC.dce = [];
  DATAC.dce_xyz_dsl = [];
  DATAC.dcv = [];
  DATAC.hk_101 = [];
  DATAC.phase = [];
  DATAC.probe2sc_pot = [];
  DATAC.sc_pot = [];
  return
end

if ~isfield(DATAC, 'scId')
  errStr = 'Data mamager not initialized! Run: mms_sdc_sdp_datamanager(''init'',init_struct)';
  irf.log('critical', errStr);
  error('MATLAB:MMS_SDC_SDP_DATAMANAGER:INPUT', errStr);
end

% If one input is give, we give the parameter back if it is already there,
% otherwise we compute it
if nargin==1
  if ~isfield(DATAC,param)
    errStr = ['unknown parameter (' param ')'];
    irf.log('critical',errStr), error(errStr)
  end
  
  if ~isempty(DATAC.(param)), varOut = DATAC.(param); return, end
  
  % Check for function fo compute param
  funcName = ['mms_sdc_sdp_comp_' param];
  if ~( exist(funcName)==2 )
    irf.log('warning',...
      ['Cannot find function (' funcName ') to compute param(' param ')'])
    DATAC.(param) = MMS_CONST.Error; varOut = DATAC.(param); return
  end
  % Run the function
  DATAC.(param) = feval(funcName); varOut = DATAC.(param); return
end

%% nargin==2   
% Make sure first argument is a dataobj class object,
% otherwise a read cdf file.
if isa(dataObj,'dataobj') % do nothing
elseif ischar(dataObj) && exist(dataObj, 'file')
  % If it is not a read cdf file, is it an unread cdf file? Read it.
  irf.log('warning',['Loading ' param ' from file: ', dataObj]);
  dataObj = dataobj(dataObj, 'KeepTT2000');
else
  errStr = 'Unrecognized input argument';
  irf.log('critical', errStr);
  error('MATLAB:MMS_SDC_SDP_DATAMANAGER:INPUT', errStr);
end

if( isfield(DATAC, param) ) && ~isempty(DATAC.(param))
  % Error, Warning or Notice for replacing the data variable?
  errStr = ['replacing existing variable (' param ') with new data'];
  irf.log('critical', errStr);
  error('MATLAB:MMS_SDC_SDP_DATAMANAGER:INPUT', errStr);
end

varPrefix = sprintf('mms%d_sdp_',DATAC.scId);

switch(param)
  case('dce')
    sensors = {'e12','e34','e56'};
    init_param()
    
  case('dcv')
    sensors = {'v1','v2','v3','v4','v5','v6'};
    init_param()
    v_from_e_and_v()
    chk_latched_p()
    chk_bias_guard()
    chk_sweep_on()
    chk_sdp_v_vals()
    
  case('hk_101')
    varPrefix = sprintf('mms%d_101_',DATAC.scId);
    DATAC.(param) = [];
    DATAC.(param).dataObj = dataObj;
    x = getdep(dataObj,[varPrefix 'cmdexec']);
    DATAC.(param).time = x.DEPEND_O.data;
    check_monoton_timeincrease(DATAC.(param).time, param);
    % Add sunpulse times (TT2000) of last recieved sunpulse.
    DATAC.(param).sunpulse = dataObj.data.([varPrefix 'sunpulse']).data;
    % Add sunpulse indicator, real: 0, SC pseudo: 1, CIDP pseudo: 2.
    DATAC.(param).sunssps = dataObj.data.([varPrefix 'sunssps']).data;
    % Add CIDP sun period (in microseconds, 0 if sun pulse not real.
    DATAC.(param).iifsunper = dataObj.data.([varPrefix 'iifsunper']).data;
  otherwise
    % Not yet implemented.
    errStr = [' unknown parameter (' param ')'];
    irf.log('critical',errStr);
    error('MATLAB:MMS_SDC_SDP_DATAMANAGER:INPUT', errStr);
end

  function chk_latched_p()
    % Check that probe values are varying. If there are 3 identical points,
    % or more, after each other mark this as latched data. If it is latched
    % and the data has a value below MMS_CONST.Limit.LOW_DENSITY_SATURATION
    % it will be Bitmasked with Low density saturation otherwise it will be
    % bitmasked with just Probe saturation.

    % Bits used for Saturation
    Bits=[MMS_CONST.Bitmask.PROBE_SATURATION, MMS_CONST.Bitmask.LOW_DENSITY_SATURATION];
    
    % For each sensor, check each pair, i.e. V_1 & V_2 and E_12.
    for iSen = 1:2:numel(sensors)
      senA = sensors{iSen};  senB = sensors{iSen+1};
      senE = ['e' senA(2) senB(2)]; % E-field sensor

      % Check all probes
      irf.log('notice', ...
        sprintf('Checking for latched probes on %s, %s and %s.', senA, ...
        senB, senE));

      indA = irf_latched_idx(DATAC.dcv.(senA).data);
      indB = irf_latched_idx(DATAC.dcv.(senB).data);
      indE = irf_latched_idx(DATAC.dce.(senE).data);

      % Add appropriate value to bitmask, leaving other 16 bits untouched.
      DATAC.dcv.(senA).bitmask(indA) = bitand(DATAC.dcv.(senA).bitmask(indA), hex2dec('FFFF')-sum(Bits)) + ...
        bitand(latched_mask(DATAC.dcv.(senA).data(indA)), sum(Bits));
      DATAC.dcv.(senB).bitmask(indB) = bitand(DATAC.dcv.(senB).bitmask(indB), hex2dec('FFFF')-sum(Bits)) + ...
        bitand(latched_mask(DATAC.dcv.(senB).data(indB)), sum(Bits));
      DATAC.dce.(senE).bitmask(indE) = bitand(DATAC.dce.(senE).bitmask(indE), hex2dec('FFFF')-sum(Bits)) + ...
        bitand(latched_mask(DATAC.dce.(senE).data(indE)), sum(Bits));

      %% TODO, Check overlapping stuck values, if senA stuck but not senB..
      
    end
    
      function latchBitmask = latched_mask( data )
        % Return value to add to bitmask, for latched probe data.
        if(~isempty(data))
          irf.log('notice', 'Latched data found');
          latchBitmask = Bits(1)*ones(size(data),'uint16');
          % Check if any of these are data with values below limit, then
          % use different latched bitmask.
          latchBitmask(data<MMS_CONST.Limit.LOW_DENSITY_SATURATION) = Bits(2);
        end
      end
  end

  function chk_bias_guard()
    % Check that bias/guard setting are nominal.
    % If not, set bit in both V and E bitmask
    
    %XXX: Does nothing at the moment
  end

  function chk_sweep_on()
    % Check if sweep is on for all probes
    % if yes, set bit in both V and E bitmask

    %XXX: Does nothing at the moment

    % Notes: cdf files (dcv/dce) will contain mmsX_sdp_sweepstatus (uint8)
    % and corresponding timestamps in mmsX_sdp_epoch_sweep (tt2000)
    % indicating if sweep is ongoing at that particular epoch time or not.
  end

  function chk_sdp_v_vals()
    % check if probe-to-spacecraft potentials  averaged over one spin for 
    % all probes are similar (within TBD %, or V). 
    % If not, set bit in both V and E bitmask.
    
    %XXX: Does nothing at the moment
  end

  function v_from_e_and_v
    % Compute V from E and the other V
    % typical situation is V2 off, V1 on
    % E12[mV/m] = ( V1[V] - V2[V] ) / L[km]
    if isempty(DATAC.dce),
      irf.log('warning','Empty DCE, cannot proceed')
      return
    end
    
    % Nominal boom length used in L1b processor
    NOM_BOOM_L = .001; % 1m, XXX the tru value should be 120 m
    
    MSK_OFF = MMS_CONST.Bitmask.SIGNAL_OFF;
    for iSen = 1:2:numel(sensors)
      senA = sensors{iSen}; senB = sensors{iSen+1};
      senE = ['e' senA(2) senB(2)]; % E-field sensor
      senA_off = bitand(DATAC.dcv.(senA).bitmask, MSK_OFF);
      senB_off = bitand(DATAC.dcv.(senB).bitmask, MSK_OFF);
      idxOneSig = xor(senA_off,senB_off);
      if ~any(idxOneSig), return, end
      iVA = idxOneSig & ~senA_off;
      if any(iVA),
        irf.log('notice',...
          sprintf('Computing %s from %s and %s for %d data points',...
          senB,senA,senE,sum(iVA)))
        DATAC.dcv.(senB).data(iVA) = DATAC.dcv.(senA).data(iVA) - ...
          NOM_BOOM_L*DATAC.dce.(senE).data(iVA);
      end
      iVB = idxOneSig & ~senB_off;
      if any(iVB),
        irf.log('notice',...
          sprintf('Computing %s from %s and %s for %d data points',...
          senA,senB,senE,sum(iVA)))
        DATAC.dcv.(senA).data(iVB) = DATAC.dcv.(senB).data(iVB) + ...
          NOM_BOOM_L*DATAC.dce.(senE).data(iVB);
      end
    end
  end

  function init_param
    DATAC.(param) = [];
    if ~all(diff(dataObj.data.([varPrefix 'samplerate_' param]).data)==0)
      err_str = 'MMS_SDC_SDP_DATAMANAGER changing sampling rate not yet implemented.';
      irf.log('warning', err_str);
      %error('MATLAB:MMS_SDC_SDP_DATAMANAGER:INPUT', err_str);
    end
    DATAC.(param).dataObj = dataObj;
    fileVersion = DATAC.(param).dataObj.GlobalAttributes.Data_version{:};
    % Skip the intial "v" and split it into [major minor revision].
    fileVersion = str2double(strsplit(fileVersion(2:end),'.'));
    DATAC.(param).fileVersion = struct('major', fileVersion(1), 'minor',...
        fileVersion(2), 'revision', fileVersion(3));
    % Make sure it is not too old to work properly.
    if DATAC.(param).fileVersion.major < MMS_CONST.MinFileVer
      err_str = sprintf('File too old: major version %d < %d',...
        DATAC.(param).fileVersion.major, MMS_CONST.MinFileVer);
      irf.log('critical',err_str), error(err_str);
    end
    x = getdep(dataObj,[varPrefix param '_sensor']);
    DATAC.(param).time = x.DEPEND_O.data;
    check_monoton_timeincrease(DATAC.(param).time, param);
    sensorData = dataObj.data.([varPrefix param '_sensor']).data;
    if isempty(sensors), return, end
    probeEnabled = resample_probe_enable(sensors);
    for iSen=1:numel(sensors)
      DATAC.(param).(sensors{iSen}) = struct(...
        'data',sensorData(:,iSen), ...
        'bitmask',zeros(size(sensorData(:,iSen)),'uint16'));
      %Set disabled bit
      idxDisabled = probeEnabled(:,iSen)==0;
      DATAC.(param).(sensors{iSen}).bitmask(idxDisabled) = ...
        bitor(DATAC.(param).(sensors{iSen}).bitmask(idxDisabled), ...
        MMS_CONST.Bitmask.SIGNAL_OFF);
      DATAC.(param).(sensors{iSen}).data(idxDisabled,:) = NaN;
    end
  end

  function res = resample_probe_enable(fields)
  % resample probe_enabled data to E-field cadense
    probe = fields{1};
    flag = get_variable(dataObj,[varPrefix probe '_enable']);
    dtSampling = median(diff(flag.DEPEND_0.data));
    switch DATAC.tmMode
%      case MMS_CONST.TmMode.srvy, error('kaboom')
      case MMS_CONST.TmMode.slow, dtNominal = [20, 160]; % seconds
      case MMS_CONST.TmMode.fast, dtNominal = 5;
      case MMS_CONST.TmMode.brst, dtNominal = [0.625, 0.229 0.0763];
      otherwise
        errS = 'Unrecognized tmMode';
        irf.log('critical',errS), error(errS)
    end
    dtNominal = int64(dtNominal*1e9); % convert to ns
    
    flagOK = false;
    for i=1:numel(dtNominal)
      if dtSampling > dtNominal(i)*.95 && dtSampling < dtNominal(i)*1.05
        dtSampling = dtNominal(i); flagOK = true; break
      end
    end
    if ~flagOK
      errS = ['bad sampling for ' varPrefix probe '_enable'];
      irf.log('critical',errS), error(errS)
    end
    enabled.time = flag.DEPEND_0.data;
    nData = numel(enabled.time);
    enabled.data = zeros(nData,numel(fields));
    enabled.data(:,1) = flag.data;
    for iF=2:numel(fields)
      probe = fields{iF};
      flag = getv(dataObj,[varPrefix probe '_enable']);
      if isempty(flag)
        errS = ['cannot get ' varPrefix probe '_enable'];
        irf.log('critical',errS), error(errS)
      elseif numel(flag.data) ~= nData
        errS = ['bad size for ' varPrefix probe '_enable'];
        irf.log('critical',errS), error(errS)
      end
      enabled.data(:,iF) = flag.data;
    end
    newT = DATAC.(param).time;
    % Default to zero - probe disabled
    res = zeros(numel(newT), numel(fields));
    if all(diff(enabled.data))==0,
      ii = newT>enabled.time(1)-dtSampling & newT<=enabled.time(end);
      for iF=1:numel(fields), 
        res(ii,iF) = enabled.data(1,iF); 
      end
    else
      % TODO: implements some smart logic.
      errS = 'MMS_SDC_SDP_DATAMANAGER enabling/disabling probes not yet implemented.';
      irf.log('critical', errS); error(errS);
    end
  end
end


% Short function for verifying Time is increasing.
function check_monoton_timeincrease(time, dataType)
    
if(any(diff(time)<=0))
        err_str = ['Time is NOT increasing for the datatype ', dataType];
        irf.log('critical', err_str);
        error('MATLAB:MMS_SDC_SDP_DATAMANAGER:TIME:NONMONOTON', err_str);
end

end
