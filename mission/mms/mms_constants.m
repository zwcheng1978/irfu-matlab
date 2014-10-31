function MMS_CONST = mms_constants
%MMS_CONST  initialize MMS constants
%
%  MMS_CONST = mms_constants()

% Version numbering, start with X, Y, Z = 0, 0, 0. When releasing new
% software update values here and subsequent output files created will have
% these numbers. 
% When simply re-running a dataset, the Z value should be increased by one.

MMS_CONST.Version.X = 1; % Major new Software version
MMS_CONST.Version.Y = 0; % New Calibration version
MMS_CONST.Version.Z = 0; % File revision, increased by 1 for each re-run.
% Version.MODS - MODS cdf GlobalAttribute should contain a description of
% all significant changes to the data set, essentially capturing a log of
% high-level release notes. Can have as many entries as necessary and
% should be updated if the "X" value of the version number changes.
% Each cell corresponds to one version, append like: mods=[mods; {'new text'}];
MMS_CONST.Version.MODS = {'V.0. Initial interface tests.'};
MMS_CONST.Version.MODS = [MMS_CONST.Version.MODS; {'V.1. Updated output variable names, added some processing.'}];


MMS_CONST.MinFileVer = 2; % min version of l1b files accepted 

MMS_CONST.MMSids = 1:4;

% Spin rate max and min, nominally 3.0 rpm +/-0.2 rpm.
MMS_CONST.Spinrate.max = 3.2; % Rev per Minute.
MMS_CONST.Spinrate.min = 2.8; % Rev per Minute.

% Telemetry mode
MMS_CONST.TmModes = {'slow','fast','brst'};
MMS_CONST.TmMode.slow = 1; % Number must corrspond to position in the list
MMS_CONST.TmMode.fast = 2;
MMS_CONST.TmMode.brst = 3;

% SDC process names
MMS_CONST.SDCProcs = {'sitl','ql','scpot','l2pre'};
MMS_CONST.SDCProc.sitl = 1; % Number must corrspond to position in the list
MMS_CONST.SDCProc.ql   = 2;
MMS_CONST.SDCProc.scpot  = 3;
MMS_CONST.SDCProc.l2pre  = 4;

% Limits used in processing
MMS_CONST.Limit.LOW_DENSITY_SATURATION = -100; % Probe stuck and below limit.

% Bitmask values; 2^(bit_number - 1):
MMS_CONST.Bitmask.SIGNAL_OFF               =  1;       % Bit 1
MMS_CONST.Bitmask.BAD_BIAS                 =  2;       % Bit 2
MMS_CONST.Bitmask.PROBE_SATURATION         =  4;       % Bit 3
MMS_CONST.Bitmask.LOW_DENSITY_SATURATION   =  8;       % Bit 4
MMS_CONST.Bitmask.SWEEP_DATA               =  16;      % Bit 5

MMS_CONST.Error = -Inf; % Indicates error in computation

% % DC V source bitmasks
% %for 16 ks/s channels, up to 6 channels at the same time:
% MMS_CONST.Source.SCM1 = 1;      % Bit 0x00 = SCM1  enable/disable
% MMS_CONST.Source.SCM1 = 2;      % Bit 0x01 = SCM2  enable/disable
% MMS_CONST.Source.SCM3 = 4;      % Bit 0x02 = SCM3  enable/disable
% MMS_CONST.Source.V1 = 8;        % Bit 0x03 = V1    enable/disable
% MMS_CONST.Source.V2 = 16;       % Bit 0x04 = V2    enable/disable
% MMS_CONST.Source.V3 = 32;       % Bit 0x05 = V3    enable/disable
% MMS_CONST.Source.V4 = 64;       % Bit 0x06 = V4    enable/disable
% MMS_CONST.Source.V5 = 128;      % Bit 0x07 = V5    enable/disable
% MMS_CONST.Source.V6 = 256;      % Bit 0x08 = V6    enable/disable
% MMS_CONST.Source.E12DC = 512;   % Bit 0x09 = E12DC enable/disable
% MMS_CONST.Source.E34DC = 1024;  % Bit 0x10 = E34DC enable/disable
% MMS_CONST.Source.E56DC = 2048;  % Bit 0x11 = E56DC enable/disable
% 
% % DC E source bitmasks
% %for 256 ks/s channels (ACE and High Speed Burst), up to 3 channels at 
% %the same time:
% MMS_CONST.Source.E12_AC = 1;    % Bit 0x00 = E12_AC enable/disable
% MMS_CONST.Source.E34_AC = 2;    % Bit 0x01 = E34_AC enable/disable
% MMS_CONST.Source.E56_AC = 4;    % Bit 0x02 = E56_AC enable/disable
% MMS_CONST.Source.V1_AC = 8;     % Bit 0x03 = V1_AC  enable/disable
% MMS_CONST.Source.V2_AC = 16:    % Bit 0x04 = V2_AC  enable/disable
