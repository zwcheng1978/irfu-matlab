classdef (Abstract) GenericTimeArray
%GenericTimeArray Generic (Abstract) class describing a time array
%
%  Methods:
%     disp()
%     end()
%     isempty()
%     length()
%     start() first point of the time array
%     stop() last point of the time array
%     tlim() Returns index and records within specified time interval
%     toEpochUnix()
%     toEpochTT2000()
  
% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------
  
  properties 
    epoch
  end
    
  methods
    function disp(obj)
      %DISP display object
      %
      %  disp(obj)
      n = obj.length();
      fprintf('  %s object : %d records\n\n',class(obj),n);
      if n>5
        for i = 1:2, dispI(i), end
        fprintf('  ... skipped %d records ...\n',n-4);
        for i = n-1:n, dispI(i), end
      else
        for i = 1:n, dispI(i), end
      end
      fprintf('\n')
      function dispI(i)
        S.type = '()'; S.subs={i};
        tmp = subsref(obj,S);
        fprintf('  %s\n',tmp.toUtc)
      end
    end
    
    function res = isempty(obj)
      %ISEMPTY True for empty array.
      res = isempty(obj.epoch);
    end
    
    function res = length(obj)
      %LENGTH number of records in array.
      res = length(obj.epoch);
    end
    
    function res = end(obj,k,~)
      %END last index in array
      switch k
        case 1
          res = length(obj);
        otherwise
          res = 1;
      end
    end
    
    function res = start(obj)
      %START  first point of the time array
      if isempty(obj)
        error('irf:GenericTimeArray:tlim:badInputs',...
            'empty input')
      end
      S.type = '()'; S.subs={1};
      res = subsref(obj,S);
    end
    
    function res = stop(obj)
      %STOP  last point of the time array
      if isempty(obj)
        error('irf:GenericTimeArray:tlim:badInputs',...
          'empty input')
      end
      S.type = '()'; S.subs={length(obj)};
      res = subsref(obj,S);
    end
    
    function [varargout] = subsref(obj,idx)
      %SUBSREF handle indexing
        switch idx(1).type
          % Use the built-in subsref for dot notation
          case '.'
            [varargout{1:nargout}] = builtin('subsref',obj,idx);
          case '()'
            tmpEpoch = builtin('subsref',obj.epoch,idx);
            [varargout{1:nargout}] = feval(class(obj),tmpEpoch);
            % No support for indexing using '{}'
          case '{}'
            error('irf:GenericTimeArray:subsref',...
              'Not a supported subscripted reference')
        end
    end
    
    function [idxLim,res] = tlim(obj,inp,mode)
      %TLIM   Returns index and records within specified time interval
      %
      % [IDX,Y] = tlim(X,LIM,[MODE])
      %
      % Where MODE can be: 
      %        'and', 0 (default)
      %        'xor', 1
      %
      %	IDX contains indexes of rows that were returned
      %
      % Y is part of the X that is within interval 
      % LIM.START <= X(:,1) < LIM.STOP for "AND" mode
      %
      % Y is part of X outside the interval for "XOR" mode:
      % X(:,1) < LIM.START & X(:,1) > LIM.STOP
      
      if isempty(inp)
        error('irf:GenericTimeArray:tlim:badInputs',...
            'empty limiting array')
      end
      if nargin<3, mode = 0; end
      if ischar(mode)
        switch lower(mode)
          case 'and'
            mode = 0;
          case 'xor'
            mode = 1;
          otherwise
            error('irf:GenericTimeArray:tlim:badInputs',...
              'MODE can be ''and'' (default) of ''xor''')
        end
      elseif isnumeric(mode)
        if ~any(mode==[0 1])
          error('irf:GenericTimeArray:tlim:badInputs',...
            'MODE can be 0 (''and'', default) of 1 (''xor''')
        end
      else
        error('irf:GenericTimeArray:tlim:badInputs',...
          'MODE can be 0 (''and'', default) of 1 (''xor''')
      end
      className = class(obj);
      lim = inp.(['to' className]);
      if nargout>1, [idxLim,res] = tlimPrivate(obj,lim,mode);
      else idxLim = tlimPrivate(obj,lim,mode);
      end
    end
    
    function res = toEpochUnix(obj)
      if isa(obj,'EpochUnix'), res = obj;
      else res = EpochUnix(toUtc(obj));
      end
    end
    function res = toEpochTT2000(obj)
      if isa(obj,'EpochTT2000'), res = obj;
      else res = EpochTT2000(toUtc(obj));
      end
    end
    
    % Anstract methods
    toUtc(obj)
    %toUtc  convert to UTC time string
  end
  
  methods (Access = private)
    function [idxLim,res] = tlimPrivate(obj,inp,mode)
      % Private version of tLim, can be reloaded
      if mode==0
        idxLim = find((obj.epoch >= inp.epoch(1)) & (obj.epoch < inp.epoch(end)));
      else
        idxLim = find((obj.epoch < inp.epoch(1)) | (obj.epoch > inp.epoch(end)));
      end
      if nargout>1,
        S.type = '()'; S.subs={idxLim};
        tmpEpoch = builtin('subsref',obj.epoch,S);
        res = feval(class(obj),tmpEpoch);
      end
    end
  end
  
  % XXX : Implement real checking
  methods (Static)
    function [ output_args ] = validate_iso_time_str( str )
      %verify_iso_time_str Summary of this function goes here
      %   Detailed explanation goes here
      
      output_args = true;
     
      if isempty(str), output_args = false; end
    end
    function res = LeapSeconds()
      % Try to read CDFLeapSeconds.txt (from NASA_cdf_patch). If not found
      % revert to hard coded values.
      % Source of both: <http://maia.usno.navy.mil/ser7/tai-utc.dat>
      % Updated: 20120401
      % Leap Seconds Table - used by CDF
      % Update it when a leap second(s) is added.
      %Year Month Day Leap Seconds      Drift
      if(exist('CDFLeapSeconds.txt', 'file'))
        n = importdata('CDFLeapSeconds.txt', ' ', 6);
        res = n.data;
        % Info about CDFLeapSeconds and time of last update.
        %info_str=['CDFLeapSeconds.txt was used.', n.textdata{2,1}(2:end)];
      else
        % File not found?
        %info_str=['No CDFLeapSeconds.txt found. Using hardcoded values,
        %updated 20120401.'];
        res = [...
          1960   1    1    1.4178180  37300.0  0.001296
          1961   1    1    1.4228180  37300.0  0.001296
          1961   8    1    1.3728180  37300.0  0.001296
          1962   1    1    1.8458580  37665.0  0.0011232
          1963  11    1    1.9458580  37665.0  0.0011232
          1964   1    1    3.2401300  38761.0  0.001296
          1964   4    1    3.3401300  38761.0  0.001296
          1964   9    1    3.4401300  38761.0  0.001296
          1965   1    1    3.5401300  38761.0  0.001296
          1965   3    1    3.6401300  38761.0  0.001296
          1965   7    1    3.7401300  38761.0  0.001296
          1965   9    1    3.8401300  38761.0  0.001296
          1966   1    1    4.3131700  39126.0  0.002592
          1968   2    1    4.2131700  39126.0  0.002592
          1972   1    1   10.0            0.0  0.0
          1972   7    1   11.0            0.0  0.0
          1973   1    1   12.0            0.0  0.0
          1974   1    1   13.0            0.0  0.0
          1975   1    1   14.0            0.0  0.0
          1976   1    1   15.0            0.0  0.0
          1977   1    1   16.0            0.0  0.0
          1978   1    1   17.0            0.0  0.0
          1979   1    1   18.0            0.0  0.0
          1980   1    1   19.0            0.0  0.0
          1981   7    1   20.0            0.0  0.0
          1982   7    1   21.0            0.0  0.0
          1983   7    1   22.0            0.0  0.0
          1985   7    1   23.0            0.0  0.0
          1988   1    1   24.0            0.0  0.0
          1990   1    1   25.0            0.0  0.0
          1991   1    1   26.0            0.0  0.0
          1992   7    1   27.0            0.0  0.0
          1993   7    1   28.0            0.0  0.0
          1994   7    1   29.0            0.0  0.0
          1996   1    1   30.0            0.0  0.0
          1997   7    1   31.0            0.0  0.0
          1999   1    1   32.0            0.0  0.0
          2006   1    1   33.0            0.0  0.0
          2009   1    1   34.0            0.0  0.0
          2012   7    1   35.0            0.0  0.0
        ];
      end
    end
  end
end