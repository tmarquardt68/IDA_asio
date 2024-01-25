function results = plot_combi_01(data,t_event,monitor_no)
% TO DO
% -
if ~exist('t_event','var')
    t_event = 0; % presentation number that should be zero on time axis
else
    t_event=t_event-1;
end
if ~exist('monitor_no','var')
    monitor_no = 3; % choose monitor
end
scrszes = get(0,'MonitorPositions');
while monitor_no > size(scrszes,1)
    monitor_no = monitor_no - 1; 
end

data_path = [fileparts(which(mfilename)) '\..\..\..\DATA\'];
data_path = 'C:\Users\Torsten\OneDrive - University College London\Documents\MATLAB\IDA_asio\DATA\'
if isstruct(data) 
    results = data;
else % parameter $data can be a filename
    try
        load([data_path data],'results')
    catch
        error('Give either $results or filename as parameter!')
    end
end

if ~isfield(results,'data') % plot combi runs the first time (monitor_combi does only save raw recordings)
    no_prmSets = size(results.stimulus.original_parameter_table,1);
    monitor_settings = monitor_init_combi_01(results);
    % load([data_path results.header.title])
    load([data_path data])
    if isfield(results,'data')
        results = rmfield(results,'data');
    end
    interval_order = results.stimulus.interval_order';
%?     interval_order = interval_order(1:end);
    for q = 1:length(who('raw*'))
        try
            eval(['wave = raw_' num2str(q) '.wave;']);
        catch % this catch is just for some older filename convention
            try
                eval(['wave = raw_' num2str(ceil(q/no_prmSets)) '_' num2str(mod(q-1,no_prmSets)+1) '.wave;']);
            catch % this catch is just for ending if raw files are missing (e.g. gp016)
                break
            end
        end
        results.data.presentation(q) = monitor_combi_01(q==1,interval_order(q),results,wave,monitor_settings,0);
    end

    % averaging per prmSet
    result_names = fieldnames(results.data.presentation);
    for prmSet = no_prmSets:-1:1
        for q=1:length(result_names)
            eval(['avg(' num2str(prmSet) ').' result_names{q}...
                '=zeros(size(results.data.presentation(1).' result_names{q} '));']);
        end
        n=0;
        for q2=1:length(results.data.presentation)
            if interval_order(q2) == prmSet
                n=n+1;
                for q=1:length(result_names)
                    eval(['avg(prmSet).'  result_names{q} '= avg(prmSet).' ...
                        result_names{q} '+results.data.presentation(q2).' result_names{q} ';']);
                end
            end
        end


        for q=1:length(result_names)
            eval(['avg(prmSet).'  result_names{q} '= avg(prmSet).' ...
                result_names{q} '/n;']);
        end
    end
    results.data.avg = avg;
    save([data_path data],'results','-append')
end

% define X_data on abcissa and plot series
x = seconds(-t_event*6):seconds(6):seconds((length(results.data.presentation)-t_event-1)*6);
h_fig = open_main_figure(results,x,scrszes(monitor_no,:));

plot_time_series_spectral_line_cronologically(results,x,h_fig);
CB_plot_prmSet(findobj(gcf,'Tag','axes_CAP'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h_fig = open_main_figure(results,x,scrsz)
h_fig = figure(get_figure_h([scrsz(1)+round(scrsz(3)/4*3),scrsz(2)+32,...
    round(scrsz(3)/4),round(scrsz(4)-110)]));
setappdata(h_fig,'x',x)
setappdata(h_fig,'results',results)
setappdata(h_fig,'scrsz',scrsz)
setappdata(h_fig,'h_main_fig',h_fig);
setappdata(gcf,'initial_call',1)
[a,b]=fileparts(results.header.title);
[~,c]=fileparts(a);
set(h_fig,'NumberTitle','off','Name',['Chrono magn series:' c ' ' b],...
    'CloseRequestFcn',@CB_close)
uicontrol('Parent',h_fig, ...
	'BusyAction','cancel', ...
	'Callback',@CB_toggle_chrono_avg, ...
	'Position',[1 1 60 15], ...
	'String','avg', ...
	'Style','pushbutton', ...
	'Tag','pushbutton_toggle', ...
	'UIContextMenu',[]);
uicontrol('Parent',h_fig, ...
	'BusyAction','cancel', ...
	'Callback',@CB_plot_prmSet, ...
	'Position',[63 1 60 15], ...
	'String','prmSet', ...
	'Style','pushbutton', ...
	'Tag','pushbutton_plot_prmSet_spctr', ...
	'UIContextMenu',[]);
uicontrol('Parent',h_fig, ...
	'BusyAction','cancel', ...
	'Callback',@CB_plot_time_series_phase, ...
	'Position',[125 1 60 15], ...
	'String','phase', ...
	'Style','pushbutton', ...
	'Tag','pushbutton_plot_time_series_phase', ...
	'UIContextMenu',[]);
uicontrol('Parent',h_fig, ...
	'BusyAction','cancel', ...
	'Callback',@CB_plot_time_series_op, ...
	'Position',[187 1 60 15], ...
	'String','op', ...
	'Style','pushbutton', ...
	'Tag','pushbutton_plot_time_series_op', ...
	'UIContextMenu',[]);
uicontrol('Parent',h_fig, ...
	'BusyAction','cancel', ...
	'Callback',@CB_delete_all_prmSet_figs, ...
	'Position',[400 1 60 15], ...
	'String','closeFigs', ...
	'Style','pushbutton', ...
	'Tag','pushbutton_closeFigs', ...
	'UIContextMenu',[]);
% uicontrol('Parent',h_fig, ...
% 	'BusyAction','cancel', ...
% 	'Callback',@CB_plot_time_series_phase_modulation, ...
% 	'Position',[315 1 60 15], ...
% 	'String','phase modulation', ...
% 	'Style','pushbutton', ...
% 	'Tag','pushbutton_plot_time_series_phase_modulation', ...
% 	'UIContextMenu',[]);
%==========================================================================

function CB_plot_prmSet(callingObj,~)
current_fig = get(callingObj,'Parent');
results = getappdata(current_fig,'results');
scrsz = getappdata(current_fig,'scrsz');
interval_order=results.stimulus.interval_order';

h_lines = get(gca,'Children');
x_data = get(h_lines(end),'XData');

% get number from mouse click position
if getappdata(current_fig,'initial_call')
    x = [];
    setappdata(current_fig,'initial_call',0)
else
    [x, ~] = ginput(1);
    x=round(x);
    if strcmp(get(callingObj,'Tag'),'pushbutton_plot_prmSet_spctr')
        zoom on
    end
    x_lim = get(gca,'XLim'); if isduration(x_lim), x_lim = seconds(x_lim); end
    if x < x_lim(1), x=x_lim(1); end
    if x > x_lim(2), x=x_lim(2); end
end
if isduration(x_data),x=seconds(x);end
kindOfPlot = get(current_fig,'Name');
if strcmp(kindOfPlot(1:3),'Avg')
    if isempty(x)
        x=x_data(floor(length(x_data)));
    end
    [~,number] = min(abs(x_data-x));
    if number > size(interval_order,1)
        number=size(interval_order,1);
    end
    data = results.data.avg;
else
    if isempty(x)
        x=x_data(floor(length(x_data)/2));
    end
    [~,number] = min(abs(x_data-x));
    if number > length(results.data.presentation)
        number=length(results.data.presentation);
    end
    data = results.data.presentation;
end
prmSet = interval_order(number);

% delete old parameter set pointer and plot new parameter set pointer in axes
if strcmp(get(current_fig,'Tag'),'chrono_series_fig')
    h_figs = findobj('Tag','chrono_series_fig');
elseif strcmp(get(current_fig,'Tag'),'avg_series_fig')
        h_figs = findobj('Tag','avg_series_fig');
end
curr_fig_name=get(current_fig,'Name');
idx=strfind(curr_fig_name,':');
for q=1:length(h_figs)
    if contains(get(h_figs(q),'Name'), curr_fig_name(idx(1):idx(1)+17))
        delete(getappdata(h_figs(q),'prmSetPointer'));
        h_ax=getappdata(h_figs(q),'h_axes');
        for q2 = 1:length(h_ax)
            axes(h_ax(q2));
            h(q2) = line([x x],ylim,'Color',[1 1 0],'Marker','o','MarkerFaceColor',[1 1 0],'MarkerSize',4);
        end
        setappdata(h_figs(q),'prmSetPointer',h)
    end
end
plot_prmSet(results,prmSet,data,number,scrsz)
hold on
%==========================================================================

function CB_toggle_chrono_avg(callingObj,~)
current_fig_name = split(get(gcbf,'Name'));
results = getappdata(gcbf,'results');
x = getappdata(gcbf,'x');
delete(getappdata(gcbf,'h_axes'))
if strcmp(current_fig_name{1},'Avg')
    plot_time_series_spectral_line_cronologically(results,x,gcbf)
    set(gcbo,'String','avg')
    set(gcbf,'Name',['Chrono ' strjoin(current_fig_name(2:end))])
else
    plot_time_series_spectral_lines_of_prmSet_averages(results,gcbf);
    set(gcbo,'String','chrono')
    set(gcbf,'Name',['Avg ' strjoin(current_fig_name(2:end))])
end
%==========================================================================

function plot_time_series_spectral_line_cronologically(results,x,h_fig)
set(h_fig,'Tag','chrono_series_fig')
result_names=fieldnames(results.data.presentation);
for q=1:length(result_names)
    eval(['series.' result_names{q} '= [results.data.presentation(:).' result_names{q} '];']);
end
plot_time_magn_series_spectral_lines(series,x,h_fig)
%==========================================================================

function plot_time_series_spectral_lines_of_prmSet_averages(results,h_fig)
set(h_fig,'Tag','avg_series_fig')
result_names = fieldnames(results.data.presentation);
for q=1:length(result_names)
    eval(['series.' result_names{q} '= [results.data.avg(:).' result_names{q} '];']);
end
% define X_data on abcissa 
x = 1:length(series.l_unmod_f2_f1);
plot_time_magn_series_spectral_lines(series,x,h_fig)
%==========================================================================

function CB_plot_time_series_phase(callingObj,~)
kindOfPlot = get(gcbf,'Name');
results = getappdata(gcbf,'results');
scrsz = getappdata(gcbf,'scrsz');
x = getappdata(gcbf,'x');
h_fig = figure(get_figure_h([scrsz(1)+round(scrsz(3)/4*2),scrsz(2)+32,...
    round(scrsz(3)/4),round(scrsz(4)-110)]));
delete(getappdata(h_fig,'h_axes'));
setappdata(h_fig,'results',results)
setappdata(h_fig,'scrsz',scrsz)
setappdata(h_fig,'h_main_fig',get(callingObj,'Parent'));
[a,b]=fileparts(results.header.title);
[~,c]=fileparts(a);
uicontrol('Parent',h_fig, ...
	'BusyAction','cancel', ...
	'Callback',@CB_plot_prmSet, ...
	'Position',[65 1 60 15], ...
	'String','prmSet', ...
	'Style','pushbutton', ...
	'Tag','pushbutton_plot_prmSet_mod', ...
	'UIContextMenu',[]);

if strcmp(kindOfPlot(1:3),'Avg')
    set(h_fig,'NumberTitle','off','Name',['Avg phase series:' c ' ' b],'Tag','avg_series_fig')
    plot_time_series_phase_of_prmSet_averages(results,h_fig);
else
    set(h_fig,'NumberTitle','off','Name',['Chrono phase series:' c ' ' b],'Tag','chrono_series_fig')
    plot_time_series_phase_cronologically(results,x,h_fig)
end
%==========================================================================

function plot_time_series_phase_cronologically(results,x,h_fig)
set(h_fig,'Tag','chrono_series_fig')
result_names = fieldnames(results.data.presentation);
for q=1:length(result_names)
    eval(['series.' result_names{q} '= [results.data.presentation(:).' result_names{q} '];']);
end
plot_time_phase_series_spectral_lines(series,x,h_fig)
%==========================================================================

function plot_time_series_phase_of_prmSet_averages(results,h_fig)
set(h_fig,'Tag','avg_series_fig')
result_names = fieldnames(results.data.presentation);
for q=1:length(result_names)
    eval(['series.' result_names{q} '= [results.data.avg(:).' result_names{q} '];']);
end
x = 1:length(series.l_unmod_f2_f1);
plot_time_phase_series_spectral_lines(series,x,h_fig)
%==========================================================================

function CB_plot_time_series_op(callingObj,~)
kindOfPlot = get(gcbf,'Name');
results = getappdata(gcbf,'results');
scrsz = getappdata(gcbf,'scrsz');
x = getappdata(gcbf,'x');
h_fig = figure(get_figure_h([scrsz(1)+round(scrsz(3)/4*2),scrsz(2)+32,...
    round(scrsz(3)/4),round(scrsz(4)-110)])); 
delete(getappdata(h_fig,'h_axes'));
setappdata(h_fig,'results',results)
setappdata(h_fig,'scrsz',scrsz)
setappdata(h_fig,'h_main_fig',get(callingObj,'Parent'));
[a,b]=fileparts(results.header.title);
[~,c]=fileparts(a);
uicontrol('Parent',h_fig, ...
	'BusyAction','cancel', ...
	'Callback',@CB_plot_prmSet, ...
	'Position',[65 1 60 15], ...
	'String','prmSet', ...
	'Style','pushbutton', ...
	'Tag','pushbutton_plot_prmSet_mod', ...
	'UIContextMenu',[]);

if strcmp(kindOfPlot(1:3),'Avg')
    set(h_fig,'NumberTitle','off','Name',['Avg OP series:' c ' ' b],'Tag','avg_series_fig')
    plot_time_series_op_of_prmSet_averages(results,h_fig);
else
    set(h_fig,'NumberTitle','off','Name',['Chrono OP series:' c ' ' b],'Tag','chrono_series_fig')
    plot_time_series_op_cronologically(results,x,h_fig)
end
%==========================================================================

function plot_time_series_op_cronologically(results,x,h_fig)
set(h_fig,'Tag','chrono_series_fig')
result_names = fieldnames(results.data.presentation);
for q=1:length(result_names)
    eval(['series.' result_names{q} '= [results.data.presentation(:).' result_names{q} '];']);
end
plot_time_op_series(series,x,h_fig)
%==========================================================================

function plot_time_series_op_of_prmSet_averages(results,h_fig)
set(h_fig,'Tag','avg_series_fig')
result_names = fieldnames(results.data.presentation);
for q=1:length(result_names)
    eval(['series.' result_names{q} '= [results.data.avg(:).' result_names{q} '];']);
end
x = 1:length(series.l_unmod_f2_f1);
plot_time_op_series(series,x,h_fig)
%==========================================================================

function CB_plot_time_series_magn_modulation(callingObj,~)
kindOfPlot = get(gcbf,'Name');
results = getappdata(gcbf,'results');
scrsz = getappdata(gcbf,'scrsz');
x = getappdata(gcbf,'x');
h_fig = figure(get_figure_h([scrsz(1)+round(scrsz(3)/4*2),scrsz(2)+32,...
    round(scrsz(3)/4),round(scrsz(4)-110)]));
delete(getappdata(h_fig,'h_axes'));
setappdata(h_fig,'results',results)
setappdata(h_fig,'scrsz',scrsz)
setappdata(h_fig,'h_main_fig',get(callingObj,'Parent'));
[a,b]=fileparts(results.header.title);
[~,c]=fileparts(a);
uicontrol('Parent',h_fig, ...
	'BusyAction','cancel', ...
	'Callback',@CB_plot_prmSet, ...
	'Position',[65 1 60 15], ...
	'String','prmSet', ...
	'Style','pushbutton', ...
	'Tag','pushbutton_plot_prmSet_mod', ...
	'UIContextMenu',[]);

if strcmp(kindOfPlot(1:3),'Avg')
    set(h_fig,'NumberTitle','off','Name',['Avg magn modulation series:' c ' ' b],'Tag','avg_series_fig')
    plot_time_magn_series_modulation(results,results.data.avg,h_fig)
else
   set(h_fig,'NumberTitle','off','Name',['Chrono magn modulation series:' c ' ' b],'Tag','chrono_series_fig')
   plot_time_magn_series_modulation(results,results.data.presentation,h_fig)
end
%==========================================================================

function CB_plot_time_series_phase_modulation(callingObj,~)
kindOfPlot = get(gcbf,'Name');
results = getappdata(gcbf,'results');
scrsz = getappdata(gcbf,'scrsz');
x = getappdata(gcbf,'x');
h_fig = figure(get_figure_h([scrsz(1)+round(scrsz(3)/4*2),scrsz(2)+32,...
    round(scrsz(3)/4),round(scrsz(4)-110)]));
delete(getappdata(h_fig,'h_axes'));
setappdata(h_fig,'results',results)
setappdata(h_fig,'scrsz',scrsz)
setappdata('h_main_fig',get(callingObj,'Parent'));
[a,b]=fileparts(results.header.title);
[~,c]=fileparts(a);
uicontrol('Parent',h_fig, ...
	'BusyAction','cancel', ...
	'Callback',@CB_plot_prmSet, ...
	'Position',[65 1 60 15], ...
	'String','prmSet', ...
	'Style','pushbutton', ...
	'Tag','pushbutton_plot_prmSet_mod', ...
	'UIContextMenu',[]);

if strcmp(kindOfPlot(1:3),'Avg')
    set(h_fig,'NumberTitle','off','Name',['Avg phase modulation series:' c ' ' b],'Tag','avg_series_fig')
    plot_time_phase_series_modulation(results,results.data.avg,x,h_fig)
else
   set(h_fig,'NumberTitle','off','Name',['Chrono phase modulation series:' c ' ' b],'Tag','chrono_series_fig')
   plot_time_phase_series_modulation(results,results.data.presentation,x,h_fig)
end

%==========================================================================

function CB_delete_all_prmSet_figs
close(findobj('Tag','prmSet_fig'))
%==========================================================================

function CB_close(callingObj,~)
scrsz = getappdata(gcbf,'scrsz');
delete(gcbf)
close(get_figure_h([scrsz(1),scrsz(2)+round(scrsz(4)*1/4)-18,round(scrsz(3)/6),round(scrsz(4)/4-40)]))
close(get_figure_h([scrsz(1),scrsz(2)+round(scrsz(4)*2/4)-18,round(scrsz(3)/6),round(scrsz(4)/4-40)]))
close(get_figure_h([scrsz(1),scrsz(2)+round(scrsz(4)*3/4)-18,round(scrsz(3)/6),round(scrsz(4)/4-40)]))
close(get_figure_h([scrsz(1)+round(scrsz(3)/6),scrsz(2)+round(scrsz(4)*2/4)-18,round(scrsz(3)/6),round(scrsz(4)/4-40)]))
close(get_figure_h([scrsz(1)+round(scrsz(3)/6*2),scrsz(2)+round(scrsz(4)*2/4)-18,round(scrsz(3)/6),round(scrsz(4)/4-40)]))
close(get_figure_h([scrsz(1)+round(scrsz(3)/6),scrsz(2)+round(scrsz(4)*3/4)-18,round(scrsz(3)/6),round(scrsz(4)/4-40)]))
close(get_figure_h([scrsz(1)+round(scrsz(3)/6*2),scrsz(2)+round(scrsz(4)*3/4)-18,round(scrsz(3)/6),round(scrsz(4)/4-40)]))
close(get_figure_h([scrsz(1)+round(scrsz(3)/4*2),scrsz(2)+32,round(scrsz(3)/4),round(scrsz(4)-110)]))

