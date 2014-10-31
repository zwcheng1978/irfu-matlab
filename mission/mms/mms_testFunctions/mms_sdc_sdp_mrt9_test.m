%% Init
data_root='/data/mms/MRT9';
%data_root='/Users/yuri/Dropbox (IRFU)/Projects/MMS/DataProcessing/Data/MRT9';
%data_root='/home/thoni/MMS/MMS_cdf/MRT9';

if ismac,
  setenv('CDF_BASE','/Applications/cdf35_0-dist/')
  outDir = '/Users/yuri/tmp_out';
else
  outDir = data_root;
end

cd(outDir)
if ~exist('log','dir'), mkdir('log'), end
if ~exist('out','dir'), mkdir('out'), end
setenv('DROPBOX_ROOT',[outDir filesep 'out'])
setenv('DATA_PATH_ROOT',[outDir filesep 'out'])
setenv('LOG_PATH_ROOT',[outDir filesep 'log'])

modes    ={'slow', 'fast', 'brst' };
versions ={'2.0.1','2.0.1','2.0.0'};
procs={'scpot','ql','sitl','l2pre'};
dates = {'20150410', '20160101'};

%% Tests
for scId=2:-1:1
  hk_101File = sprintf('%s/mms%d/mms%d_fields_hk_l1b_101_20150410_v0.1.0.cdf',...
    data_root,scId,scId);
  for iMode=1:length(modes)
    for iDate = 1:numel(dates)
      dceFile = sprintf('%s/mms%d/mms%d_sdp_%s_l1b_dce_%s*_v%s.cdf',...
        data_root,scId,scId,modes{iMode},dates{iDate},versions{iMode});
      dcvFile = sprintf('%s/mms%d/mms%d_sdp_%s_l1b_dcv_%s*_v%s.cdf',...
        data_root,scId,scId,modes{iMode},dates{iDate},versions{iMode});
      d = dir(dceFile);
      if isempty(d), dceFile = '';
      else dceFile = sprintf('%s/mms%d/%s',data_root,scId,d(1).name);
      end
      d = dir(dcvFile);
      if isempty(d), dcvFile = '';
      else dcvFile = sprintf('%s/mms%d/%s',data_root,scId,d(1).name);
      end
      if isempty(dcvFile) && isempty(dceFile), continue, end
      for iProc=1:length(procs)
        fprintf('TEST: %s MMS%d %s %s\n',...
          dates{iDate},scId,modes{iMode},procs{iProc})
        mms_sdc_sdp_proc(procs{iProc},dceFile,dcvFile,hk_101File)
      end
    end
  end
end
