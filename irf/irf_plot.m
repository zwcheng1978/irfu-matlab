function c=irf_plot(varargin)
%IRF_PLOT   Flexible plotting routine for time series
%
% IRF_PLOT(X) plot data X
%   X is one of:
%      - matrix of data (1st column time in isdat epoch, other columns data)
%      - cell array data, each of cells containing a matrix of data 
%      - string defining variable (can be CAA variable)
%      - number (initialize figure with so many subplots)
%
% H=IRF_PLOT(..) return handle to axes of plots into C
% 
% IRF_PLOT(H,X,..) plot into axes handle H
%
% IRF_PLOT(X,arguments) arguments can be:
%     all arguments allowed by PLOT command 
%     'subplot' - plot all x values in separate subplots
%     'comp'    - plot vector component in separate subplots
%     ['dt', [dt1, dt2, dt3, dt4]] - specify time shifts, new time = old time - dt
%     ['yy',factor_to_multiply] - add second axis on right, miltiply by factor_to_multiply
%     ['linestyle',LineStyle] - define line style. Simple LineStyle can be be
%     given as last argument, 'linestyle' keyword is not necessary in this
%     case. LineStyle can be given as cell array to specify style for different variables/subplots.
%
% To improve zooming of high time resolution data IRF_PLOT will sometimes 
% set t_start_epoch within the 'userdata' field of the figure and 
% internally use it as origo but in most cases you should not care about this.
%
% IRF_PLOT(CAA_variable,..) can have all arguments as DATAOBJ/PLOT
%
% H=IRF_PLOT(number)	reset current figure to 'number' subplots
%						keep figure properties, return handle vector to subplots 
% IRF_PLOT(number,'reset') reset current figure 
%						but update default properties except figure size
% IRF_PLOT(number,'newfigure') initialize new figure with 'number' of subplots 
%
% Examples:
%   irf_plot(B1) - plot variable B1 (all components), assuming that the first column is time
%   irf_plot('B1') - plot variable B1, if it does not exist try to load it
%                    with c_load('B1') and try to put ylabel from c_desc('B1')
%   irf_plot('B?') - Cluser oriented, plot B1.. B4 in separate subplots
%   irf_plot({B1,B2}) - plot B1 and B2 in separate subplots
%   irf_plot('B1 B2') - -"- but if B1,B2 do not exist try to load them
%   irf_plot({B1,B2},'comp') - plot in 1. subplot B1_X and B2_X, in second
%                    subplot B1_Y and B2_Y etc.
%   irf_plot({B1,B2},'dt',[dt1 dt2]) - separate subplots with B1 and B2,
%                    but in addition B1 and B2 time axis are shifted by dt1
%                    and dt2 correspondingly
%   irf_plot('B_mag__C1_CP_FGM_5VPS') plot CAA variable B_mag__C1_CP_FGM_5VPS (if necessary load it)
%
% See also DATAOBJ/PLOT, C_PL_TX, C_DESC

% flag_subplot 0 - one plot
%              1 - separate subplots for every component
%              2 - separate subplots for all variables in the cell array
%              3 - components of vectors in separate panels
%
% For zooming to work even in cases of wide band it is important that time
% axis is not big number. Isdat epoch is too big. Therefore if time is
% isdat epoch we choose reference time the first point of first variable
% (in practices it does not matter).

%% Check input 
[ax,args,nargs] = axescheck(varargin{:});
x=args{1};
if isempty(x), % nothing to plot, first input parameter empty
  return;
end

% Check if single number argument, then use syntax IRF_PLOT(number)
if isnumeric(x),    
    if numel(x)==1, % only one number
		if x>=1 && x<=20,
			% check if there is 'newfigure' argument
			if numel(args)>=2 && ischar(args{2}) && strcmpi(args{2},'newfigure')
				c=initialize_figure(x,'newfigure');
			elseif numel(args)>=2 && ischar(args{2}) && strcmpi(args{2},'reset')
				c=initialize_figure(x,'reset');
			else
				c=initialize_figure(x);
			end
		else
			disp('Only 1-20 number of subplots supported.;)');
		end
		if nargout==0, clear c; end % if no output required do not return anything
        return
    end
end
if isempty(ax),
    inp_name=inputname(1);
else
    inp_name=inputname(2);
end
args=args(2:end);
original_args=args;

var_desc{1} = '';
flag_subplot = 0;
have_options = 0;
caa_dataobject={[]}; % by default assume we are not working with CAA variables

if nargs > 1, have_options = 1; end

%% Default values that can be override by options
dt = 0;
flag_yy = 0;
scaleyy = 1;
plot_type = '';
marker = '-';
flag_plot_all_data=1;
flag_colorbar=1;
%% Check input options
while have_options
    l = 1;
    switch(lower(args{1}))
        case 'newfigure'
            c=initialize_figure(x);
        case 'subplot'
            plot_type = 'subplot';
        case 'comp'
            plot_type = 'comp';
        case 'dt'
            if nargs>1
                if isnumeric(args{2})
                    dt = args{2};
                    l = 2;
                else 
					irf.log('critical','wrongArgType : dt must be numeric')
					error('dt must be numeric');
                end
            else irf.log('critical','wrongArgType : dt value is missing')
				error('dt value missing');
            end
        case 'tint'
            if nargs>1 && isnumeric(args{2})
                    tint = args{2};
                    l = 2;
                    flag_plot_all_data=0;
                    original_args(find(strcmpi(original_args,'tint'))+ [0 1])=[]; % remove tint argument
            else % TODO implement string tint
				irf.log('critical','wrongArgType : tint must be numeric')
				error('tint must be numeric')
            end
        case 'yy'
            if nargs>1
                if isnumeric(args{2})
                    flag_yy = 1;
                    scaleyy = args{2};
                    l = 2;
                else 
					irf.log('critical','wrongArgType : yy must be numeric')
					error('yy must be numeric');
                end
            else 
				irf.log('critical','wrongArgType : yy value is missing')
				error('yy value missing');
            end
        case 'linestyle'
            marker = args{2};
            l = 2;
        case 'nocolorbar'
            flag_colorbar=0;
            l = 1;
        otherwise
            marker = args{1};
            args = args(2:end);
            break
    end
    args = args(l+1:end);
    if isempty(args), break, end
end

%% Plot separate subplots for all x components
if strcmp(plot_type,'subplot') && isnumeric(x), flag_subplot = 1; end
if ischar(x), % Try to get variable labels etc.
    var_nam = tokenize(x); % White space separates variables
    jj = 1;var_names=cell(1,4*length(var_nam));
    for ii=1:length(var_nam), % construct varibale names var_names
        if regexp(var_nam{ii},'?'),
            c_eval(['var_names{jj}=''' var_nam{ii} ''';jj=jj+1;']);
        else
            var_names{jj} = var_nam{ii}; jj=jj+1;
        end
    end
    var_names(jj:end)=[];
    x = cell(1,length(var_names)); % preallocate expected number of variables
    caa_varname=x;                 % preallocate caa variable names
    var_desc=x;                    % preallocate non-caa variable description
    ix = 1;
    for ii=1:length(var_names) % get variables
        if evalin('caller',['exist(''' var_names{ii} ''')']), % Try to get variable from calling workspace
            x{ix} = evalin('caller',var_names{ii});
        elseif strfind(var_names{ii},'__') % CAA variable
                    caa_varname{ix}=var_names{ii};
                    if flag_plot_all_data,
                      [tmp,caa_dataobject{ix},x{ix}]=evalin('caller',['c_caa_var_get(''' var_names{ii} ''')']);
                    else
                      [tmp,caa_dataobject{ix},x{ix}]=c_caa_var_get(var_names{ii},'tint',tint);
                    end
        else
            flag_ok=c_load(var_names{ii});
            if flag_ok
                eval(['x{ix}=' var_names{ii} ';']);
            else
                irf.log('warning',...
                    ['skipping, do not know where to get variable >'...
                    var_names{ii}]);
            end
        end
        if ~isempty(x{ix}), % has succeeded to get new variable
            var_desc{ix} = c_desc(var_names{ii});
            ix = ix +1;
        end
    end
    x(ix:end)=[]; % remove unused variables
    caa_varname(ix:end)=[];
    var_desc(ix:end)=[];
end

if iscell(x), % Plot several variables
    
    % No ylabels are given
    % But no way to now the name of variables
    if size(var_desc,2)<size(x,2), var_desc = cell(1,length(x)); end
    
    if dt==0, dt(1:size(x,2)) = double(0); end
    
    switch plot_type
        case ''
            flag_subplot = 2;
            if length(x)==1, x = x{1}; flag_subplot = 0; end
        case 'comp'
            flag_subplot = 3;
        case 'subplot'
            flag_subplot = 2;
    end
else
    if isempty(inp_name)
        var_desc{1} = {};
    else    
        try % try to obtain variable description
            var_desc{1} = c_desc(inp_name);
        catch %#ok<CTCH>
            var_desc{1} = {};
        end
    end
end
if ~isempty(caa_dataobject{1}) % plot CAA variable
    if isempty(ax)
        ax = plot(caa_dataobject{1},caa_varname{1},original_args{:});
    else
        ax = plot(ax,caa_dataobject{1},caa_varname{1},original_args{:});
    end
    if isstruct(x), firstTimeStamp=x.t(1);
    elseif iscell(x), firstTimeStamp=x{1}(1,1);
    else firstTimeStamp=x(1,1);
    end
    c=ax; % axis to which apply irf_timeaxis
    flag_subplot=-1; % dont make more plots
end

if isempty(ax), % if empty axis use current axis GCA
    if isempty(get(0,'CurrentFigure')), % there is no figure open
        irf_plot(1);
    end
    ax=gca;
end

%% One subplot only
if flag_subplot==0,  % One subplot
    if isstruct(x)
        % Plot a spectrogram
        irf_spectrogram(ax,x);
		if flag_colorbar, hcbar = colorbar('peer',ax); end
        if ~isempty(var_desc{1})
            lab = cell(1,length(var_desc{1}.size));
            for v = 1:length(var_desc{1}.size)
                lab{v} = [var_desc{1}.labels{v} '[' var_desc{1}.units{v} ...
                    '] sc' var_desc{1}.cl_id];
            end
            ylabel(hcbar, lab);
        end
        
        firstTimeStamp = x.t(~isnan(x.t),1);
        firstTimeStamp = firstTimeStamp(1);
    elseif ~isempty(x) % x is nonempty matrix
        ts = t_start_epoch(x(:,1)); % t_start_epoch is saved in figures user_data variable
        ii = 2:length(x(1,:));
        tag=get(ax,'tag');
        ud=get(ax,'userdata'); % keep userdata during plotting
        if flag_yy == 0,
            h = plot(ax,(x(:,1)-ts-dt),x(:,ii),marker,args{:});
        else
            h = plotyy(ax,(x(:,1)-ts),x(:,ii),(x(:,1)-ts),x(:,ii).*scaleyy);
        end
        grid(ax,'on');
        set(ax,'tag',tag);
        set(ax,'userdata',ud);
        zoom_in_if_necessary(ax);
        
        % Put YLimits so that no labels are at the end (disturbing in multipanel plots)
        if ~ishold(ax), irf_zoom(ax,'y'); end % automatic zoom only if hold is not on
        
        if ~isempty(var_desc{1}) && isfield(var_desc{1},'size')
            lab = cell(1,length(var_desc{1}.size));
            for v = 1:length(var_desc{1}.size)
                lab{v} = [var_desc{1}.labels{v} '[' var_desc{1}.units{v} ...
                    '] sc' var_desc{1}.cl_id];
            end
            ylabel(ax,lab);
        end
        c=h;
%        c = get(h(1),'Parent');
        
        firstTimeStamp = x(~isnan(x(:,1)),1);
        firstTimeStamp = firstTimeStamp(1);
    else % empty matrix or does not know what to do
        return
    end
    
elseif flag_subplot==1, % Separate subplot for each component
    if isstruct(x), error('cannot plot spectra in COMP mode'), end
    
    % t_start_epoch is saved in figures user_data variable
    ts = t_start_epoch(x(:,1));
    
    npl = size(x,2) -1;
    c = zeros(1,npl);
    for ipl=1:npl
        c(ipl) = subplot(npl,1,ipl);
        
        if iscell(marker)
            if length(marker)==npl, marker_cur = marker{ipl};
            else marker_cur = marker{1};
            end
        else marker_cur = marker;
        end
        
        plot((x(:,1)-ts-dt),x(:,ipl+1),marker_cur,args{:}); grid on;
        
        % Put YLimits so that no labels are at the end (disturbing in multipanel plots)
        irf_zoom(c(ipl),'y');
        
        if ~isempty(var_desc) && ~isempty(var_desc{1})
            scu = cumsum(var_desc{1}.size);
            isz = find( scu == min(scu(ipl<=scu)) );
            sz = var_desc{1}.size(isz); % Size of a data vector
            if sz == 1 % Scalar data
                lab = [var_desc{1}.labels{isz} ' ['...
                    var_desc{1}.units{isz} '] sc' var_desc{1}.cl_id];
            else % Vector data
                % Vector component
                if isz==1, comp = ipl;
                else comp = ipl -scu(isz-1);
                end
                lab = [var_desc{1}.labels{isz} ...
                    '_{' var_desc{1}.col_labels{isz}{comp} '} ['...
                    var_desc{1}.units{isz} '] sc' var_desc{1}.cl_id ];
            end
            ylabel(lab);
        end
    end
    
    firstTimeStamp = x(~isnan(x(:,1)),1);
    firstTimeStamp = firstTimeStamp(1);
    
elseif flag_subplot==2, % Separate subplot for each variable
    if isempty(x), return, end
    
    %   t_start_epoch is saved in figures user_data variable
    if isstruct(x{1}), ts = t_start_epoch(x{1}.t);
    else ts = t_start_epoch(x{1}(:,1));
    end
    
    t_st = []; t_end = [];
    xlen = [];
    
    npl = size(x,2);
    c=initialize_figure(npl);
    for ipl=1:npl
        y = x{ipl};
        if isstruct(y), t_tmp = double(y.t);
        else t_tmp = double(y(:,1));
        end
        if numel(t_tmp)==0
            irf.log('critical',['Can not plot data ' num2str(ipl)]);
            return;
        end
        t_tmp = t_tmp -double(ts) -double(dt(ipl));
        firstTimeStamp = t_tmp(~isnan(t_tmp));
        if isempty(t_st), t_st = firstTimeStamp(1);
        else if firstTimeStamp(1)<t_st, t_st = firstTimeStamp(1); end
        end
        if isempty(t_end), t_end = firstTimeStamp(end);
        else if firstTimeStamp(end)>t_end, t_end = firstTimeStamp(end); end
        end
        clear tt
        
        if isstruct(y)
            irf_spectrogram(c(ipl),y.t-dt(ipl), y.p, y.f);
            if flag_colorbar, hcbar = colorbar; end
            if ~isempty(var_desc{ipl})
                lab = cell(1,length(var_desc{ipl}.size));
                for v = 1:length(var_desc{ipl}.size)
                    lab{v} = [var_desc{ipl}.labels{v} '[' var_desc{ipl}.units{v} ...
                        '] sc' var_desc{ipl}.cl_id];
                end
                ylabel(hcbar, lab);
                disp(lab)
            end
            firstTimeStamp = y.t(~isnan(y.t),1);
            % Save panel width to resize the rest of the panels accordingly
            if isempty(xlen)
                xlen = get(c(ipl),'Position');
                xlen = xlen(3);
            end
        else
            if iscell(marker)
                if length(marker)==npl, marker_cur = marker{ipl};
                else marker_cur = marker{1};
                end
            else marker_cur = marker;
            end
            plot(c(ipl),t_tmp,y(:,2:end),marker_cur); 
            grid(c(ipl),'on');
            zoom_in_if_necessary(c(ipl));
            
            % Put YLimits so that no labels are at the end (disturbing in multipanel plots)
            irf_zoom(c(ipl),'y');
            
            if ~isempty(var_desc) && ~isempty(var_desc{ipl})
                for v = 1:length(var_desc{ipl}.size)
                    lab{v} = [var_desc{ipl}.labels{v} '[' ...
                        var_desc{ipl}.units{v} '] sc' var_desc{ipl}.cl_id];
                end
                ylabel(c(ipl),lab); clear lab
            end
            firstTimeStamp = y(~isnan(y(:,1)),1);
        end
    end
    % Set common XLim
    for ipl=1:npl
        set(c(ipl),'XLim',[t_st t_end])
        if ~isempty(xlen)
            p = get(c(ipl),'Position');
            set(c(ipl),'Position',[p(1) p(2) xlen p(4)])
        end
    end
    clear t_st t_end
    
    firstTimeStamp = firstTimeStamp(1);
    
elseif flag_subplot==3,  % components of vectors in separate panels
    if isstruct(x), error('cannot plot spectra in COMP mode'), end
    % t_start_epoch is saved in figures user_data variable
    ts = t_start_epoch(x{1}(:,1));
    
    npl = size(x{1},2) -1;
    if npl==1,     % We make new figure with subplots only if more than 1 component to plot
        c = ax;
    else
        c=initialize_figure(npl);
    end
    for ipl=1:npl
        
        line_colors=get(c(ipl),'ColorOrder');
        for jj=1:size(x,2)
            use_color = 1;
            if iscell(marker)
                if length(marker)==size(x,2), marker_cur = marker{jj};  use_color = 0;
                else marker_cur = marker{1};
                end
            else marker_cur = marker;
            end
            
            if size(x{jj},2)>=ipl+1
                y = x{jj};
                if use_color
                    plot(c(ipl),(y(:,1)-ts-dt(jj)), y(:,ipl+1),...
                        'Color', line_colors(jj,:), 'LineStyle',marker_cur)
                else
                    plot(c(ipl),(y(:,1)-ts-dt(jj)), y(:,ipl+1),marker_cur)
                end
                hold(c(ipl),'on');
            end
        end
        grid(c(ipl),'on');
        
        % Put YLimits so that no labels are at the end (disturbing in
        % multipanel plots)
        set(c(ipl),'YLim',...
            mean(get(c(ipl),'YLim'))+diff(get(c(ipl),'YLim'))*[-.499999 .499999])
        
    end
    firstTimeStamp = y(~isnan(y(:,1)),1);
    firstTimeStamp = firstTimeStamp(1);
end

%% Add figure menu
irf_figmenu;

%% Add information about subplot handles to userdata of figure
%user_data = get(gcf,'userdata');
%if flag_subplot>0, user_data.subplot_handles = c; end
%set(gcf,'userdata',user_data);

%% In case time is in isdat_epoch add time axis
if ((firstTimeStamp > 1e8) && (firstTimeStamp < 1e10))
    if flag_subplot == 0, irf_timeaxis(ax);
	elseif strcmp(get(c(1),'type'),'axes')
		irf_timeaxis(c);
	elseif strcmp(get(get(c(1),'parent'),'type'),'axes')
		irf_timeaxis(get(c(1),'parent'));
	else % do nothing
    end
end

%% Do not give axis handle as answer if not asked for
if nargout==0, clear c; end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function t_st_e = t_start_epoch(t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gives back the value of t_start_epoch of the figure
% if not  set, sets t_start_epoch of the figure
ud = get(gcf,'userdata');
ii = find(~isnan(t));
if ~isempty(ii), valid_time_stamp = t(ii(1)); else valid_time_stamp = []; end

if isfield(ud,'t_start_epoch')
    t_st_e = double(ud.t_start_epoch);
elseif ~isempty(valid_time_stamp)
    if valid_time_stamp > 1e8
        % Set start_epoch if time is in isdat epoch
        % Warn about changing t_start_epoch
        t_st_e = double(valid_time_stamp);
        ud.t_start_epoch = t_st_e;
        set(gcf,'userdata',ud);
        irf.log('notice',['user_data.t_start_epoch is set to ' ...
            epoch2iso(t_st_e,1)]);
    else
        t_st_e = double(0);
    end
else
    t_st_e = double(0);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c=initialize_figure(number_of_subplots,flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% flag = "newfigure" % [optional] if to open a new figure
if nargin==1, flag='';end
if isempty(get(0,'CurrentFigure')) % no current figures opened
    flag='newfigure';
elseif isempty(get(gcf,'children')) && ~strcmpi(flag,'newfigure') && ~strcmpi(flag,'reset') % current figure is empty, display warning of new syntax (REMOVE THIS ELSE LOOP IN 2013)
	disp('WARNING! use syntax irf_plot(number_of_subplots,''newfigure'') if you want new figure.')
end
if number_of_subplots>=1 && number_of_subplots<=20,
    number_of_subplots=floor(number_of_subplots);
    c=zeros(1,number_of_subplots);
	if strcmpi(flag,'newfigure'), % if to open new figure
		hcf = figure;
		xSize = 11;
		ySize = 5+5*sqrt(number_of_subplots);
		xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
		set(hcf,'PaperPosition',[xLeft yTop xSize ySize])
		un=get(0,'units');
		set(0,'units','pixels');
		sz=get(0,'screensize');
		xx=min(min(700,sz(3))/xSize,min(900,sz(4))/ySize); % figure at least 600 wide or 900 height but not outside screen
		set(hcf,'Position',[10 10 xSize*xx ySize*xx])
		set(0,'units',un);
		clear xSize sLeft ySize yTop
	else
		hcf = gcf; 
	end
	if strcmpi(flag,'newfigure') || strcmpi(flag,'reset') 
		set(hcf,'color','white'); % white background for figures (default is grey)
        set(hcf,'renderer','zbuffer'); % opengl has problems on Mac (no log scale in spectrograms)
        set(hcf,'PaperUnits','centimeters');
        set(hcf,'defaultlinelinewidth',1.0);
        set(hcf,'defaultAxesFontSize',14);
        set(hcf,'defaultTextFontSize',14);
        set(hcf,'defaultAxesFontUnits','pixels');
        set(hcf,'defaultTextFontUnits','pixels');
        set(hcf,'defaultAxesColorOrder',[0 0 0;0 0 1;1 0 0;0.3 0.3 0.3;0 1 1 ;1 0 1; 1 1 0])
	end
    clf;
    all_axis_position=[0.17 0.1 0.9 0.95]; % xmin ymin xmax ymax
    subplot_width=all_axis_position(3)-all_axis_position(1);
    subplot_height=(all_axis_position(4)-all_axis_position(2))/number_of_subplots;
    for j=1:number_of_subplots,
        c(j)=axes('position',[all_axis_position(1) ...
            all_axis_position(4)-j*subplot_height ...
            subplot_width subplot_height]); % [x y dx dy]
%        c(j)=irf_subplot(number_of_subplots,1,-j);
        cla(c(j));
        set(c(j),'box','on');
        set(c(j),'tag','');
    end
    user_data = get(gcf,'userdata');
    user_data.subplot_handles = c;
    user_data.current_panel=0;
    set(hcf,'userdata',user_data);
    figure(hcf); % bring figure to front
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function zoom_in_if_necessary(h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ud=get(h,'userdata');
if isfield(ud,'zoom_x'),
    irf.log('debug','zooming in the updated plot')
    irf_zoom(h,'x',ud.zoom_x);
    if ud.zoom_x(1) > 1e8 && ud.zoom_x(1) < 1e10, % isdat epoch
        irf_timeaxis(h,'nolabel');
    end
end
end
