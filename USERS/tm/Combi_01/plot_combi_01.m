function results = plot_combi_01(data)
% TO DO
% -

monitor_no = 3; % choose monitor

scrsz = get(0,'MonitorPositions');
while monitor_no > size(scrsz,1)
    monitor_no = monitor_no - 1; 
end

% data_path = [fileparts(which(mfilename)) '\..\..\..\DATA\'];
data_path = 'C:\Users\Admin\Documents\MATLAB\IDA_asio_June23\DATA\';
if isstruct(data) %resutls can eb a filename
    results = data;
else
    try
        load([data_path data],'results')
    catch
        error('Give either $results or filename as parameter!')
    end
end

%if ~isfield(results,'data') % plot combi runs the first time (monitor_combi does only save raw recordings)
    no_prmSets = size(results.stimulus.original_parameter_table,1);
    monitor_settings = monitor_init_combi_01(results);
    % load([data_path results.header.title])
    load([data_path data])
    if isfield(results,'data')
        results = rmfield(results,'data');
    end
    interval_order = results.stimulus.interval_order';
    interval_order = interval_order(1:end);
    for q = 1:length(who('raw*'))
        try
            eval(['wave = raw_' num2str(q) '.wave;']);
        catch % this catch is just for spme older filename convention
            eval(['wave = raw_' num2str(ceil(q/no_prmSets)) '_' num2str(mod(q-1,no_prmSets)+1) '.wave;']);
        end
        results.data.presentation(q) = monitor_combi_01(q==1,interval_order(q),results,wave,monitor_settings,0);
    end
    save([data_path data],'results','-append')
%end
result_names = fieldnames(results.data.presentation);
no_prmSets = size(results.stimulus.original_parameter_table,1);
interval_order = results.stimulus.interval_order';
interval_order = interval_order(1:end);

% !!! THIS CAN GO LATER INSIDE IF BRANCH ABOVE
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


plot_time_series_spectral_line_cronologically(results,scrsz(monitor_no,:));
setappdata(gcf,'initial_call',1)
CB_plot_prmSet

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CB_plot_prmSet(callingObj,~)
results = getappdata(gcf,'results');
scrsz = getappdata(gcf,'scrsz');

% get number from mouse click position
if getappdata(gcf,'initial_call')
    number = [];
    setappdata(gcf,'initial_call',0)
else
    [x, ~] = ginput(1);
    zoom on
    number = round(x);
    if number < 1
        number=1;
    end
    if number > size(results.stimulus.original_parameter_table,1)
        number=size(results.stimulus.original_parameter_table,1);
    end

end

kindOfPlot = get(gcf,'Name');
if strcmp(kindOfPlot(1:3),'Avg')
    if isempty(number)
        number = round(size(results.stimulus.original_parameter_table,1)/2);
    end
    plot_prmSet(results,number,results.data.avg,number,scrsz)
else
    interval_order = results.stimulus.interval_order';
    if isempty(number)
        number = round(length(results.data.presentation)/2);
    end
    prmSet = interval_order(number);
    plot_prmSet(results,prmSet,results.data.presentation,number,scrsz)
end
%==========================================================================

function CB_toggle_chrono_avg(callingObj,~)
kindOfPlot = get(gcbf,'Name');
results = getappdata(gcbf,'results');
scrsz = getappdata(gcbf,'scrsz');

if strcmp(kindOfPlot(1:3),'Avg')
    plot_time_series_spectral_line_cronologically(results,scrsz)
    set(findobj('Tag','pushbutton_toggle'),'String','plot_avg')
else
    plot_time_series_spectral_lines_of_prmSet_averages(results,scrsz);
    set(findobj('Tag','pushbutton_toggle'),'String','plot_chrono')
end
%==========================================================================

function plot_time_series_spectral_line_cronologically(results,scrsz)
h_fig = figure(get_figure_h([scrsz(1)+round(scrsz(3)/4*3),scrsz(2)+52,...
    round(scrsz(3)/4),round(scrsz(4)-110)])); clf
setappdata(h_fig,'results',results)
setappdata(h_fig,'scrsz',scrsz)
set(h_fig,'Name',['Chrono series' results.header.title])
uicontrol('Parent',h_fig, ...
	'BusyAction','cancel', ...
	'Callback',@CB_toggle_chrono_avg, ...
	'Position',[1 1 60 15], ...
	'String','plot avg', ...
	'Style','pushbutton', ...
	'Tag','pushbutton_toggle_chrono_avg', ...
	'UIContextMenu',[]);
uicontrol('Parent',h_fig, ...
	'BusyAction','cancel', ...
	'Callback',@CB_plot_prmSet, ...
	'Position',[63 1 60 15], ...
	'String','plot prmSet', ...
	'Style','pushbutton', ...
	'Tag','pushbutton_plot_prmSet', ...
	'UIContextMenu',[]);
uicontrol('Parent',h_fig, ...
	'BusyAction','cancel', ...
	'Callback',@CB_plot_time_series_phase, ...
	'Position',[125 1 60 15], ...
	'String','plot phase', ...
	'Style','pushbutton', ...
	'Tag','pushbutton_plot_time_series_phase', ...
	'UIContextMenu',[]);
uicontrol('Parent',h_fig, ...
	'BusyAction','cancel', ...
	'Callback',@CB_plot_time_series_modulation, ...
	'Position',[315 1 80 15], ...
	'String','plot modulation', ...
	'Style','pushbutton', ...
	'Tag','pushbutton_plot_time_series_modulation', ...
	'UIContextMenu',[]);

result_names=fieldnames(results.data.presentation);
for q=1:length(result_names)
    eval(['series.' result_names{q} '= [results.data.presentation(:).' result_names{q} '];']);
end
plot_time_magn_series_spectral_lines(series)
%==========================================================================

function plot_time_series_spectral_lines_of_prmSet_averages(results,scrsz)
h_fig = figure(get_figure_h([scrsz(1)+round(scrsz(3)/4*3),scrsz(2)+52,...
    round(scrsz(3)/4),round(scrsz(4)-110)])); clf
setappdata(h_fig,'results',results)
setappdata(h_fig,'scrsz',scrsz)
set(h_fig,'Name',['Avg series' results.header.title])
uicontrol('Parent',h_fig, ...
	'BusyAction','cancel', ...
	'Callback',@CB_toggle_chrono_avg, ...
	'Position',[1 1 60 15], ...
	'String','plot chrono', ...
	'Style','pushbutton', ...
	'Tag','pushbutton_toggle', ...
	'UIContextMenu',[]);
uicontrol('Parent',h_fig, ...
	'BusyAction','cancel', ...
	'Callback',@CB_plot_prmSet, ...
	'Position',[63 1 60 15], ...
	'String','plot prmSet', ...
	'Style','pushbutton', ...
	'Tag','pushbutton_plot_prmSet', ...
	'UIContextMenu',[]);
uicontrol('Parent',h_fig, ...
	'BusyAction','cancel', ...
	'Callback',@CB_plot_time_series_phase, ...
	'Position',[125 1 60 15], ...
	'String','plot phase', ...
	'Style','pushbutton', ...
	'Tag','pushbutton_plot_time_series_phase', ...
	'UIContextMenu',[]);
uicontrol('Parent',h_fig, ...
	'BusyAction','cancel', ...
	'Callback',@CB_plot_time_series_modulation, ...
	'Position',[315 1 80 15], ...
	'String','plot modulation', ...
	'Style','pushbutton', ...
	'Tag','pushbutton_plot_time_series_modulation', ...
	'UIContextMenu',[]);

result_names = fieldnames(results.data.presentation);
for q=1:length(result_names)
    eval(['series.' result_names{q} '= [results.data.avg(:).' result_names{q} '];']);
end
plot_time_magn_series_spectral_lines(series)
%==========================================================================

function CB_plot_time_series_modulation(callingObj,~)
kindOfPlot = get(gcbf,'Name');
results = getappdata(gcbf,'results');
scrsz = getappdata(gcbf,'scrsz');

if strcmp(kindOfPlot(1:3),'Avg')
    plot_time_series_modulation_of_prmSet_averages(results,scrsz);
else
    plot_time_series_modulation_cronologically(results,scrsz)
end
%==========================================================================

function plot_time_series_modulation_cronologically(results,scrsz)
h_fig = figure(get_figure_h([scrsz(1)+round(scrsz(3)/4*2),scrsz(2)+52,...
    round(scrsz(3)/4),round(scrsz(4)-110)])); clf
setappdata(h_fig,'results',results)
setappdata(h_fig,'scrsz',scrsz)
set(h_fig,'Name',['Chrono modulation series' results.header.title])
uicontrol('Parent',h_fig, ...
	'BusyAction','cancel', ...
	'Callback',@CB_plot_prmSet, ...
	'Position',[65 1 60 15], ...
	'String','plot prmSet', ...
	'Style','pushbutton', ...
	'Tag','pushbutton_plot_prmSet', ...
	'UIContextMenu',[]);
plot_time_series_modulation(results.data.presentation)
%==========================================================================

function plot_time_series_modulation_of_prmSet_averages(results,scrsz)

t = 'Function plot_time_series_modulation_of_prmSet_averages not implemented yet'

%==========================================================================

function CB_plot_time_series_phase(callingObj,~)
kindOfPlot = get(gcbf,'Name');
results = getappdata(gcbf,'results');
scrsz = getappdata(gcbf,'scrsz');

if strcmp(kindOfPlot(1:3),'Avg')
    plot_time_series_phase_of_prmSet_averages(results,scrsz);
else
    plot_time_series_phase_cronologically(results,scrsz)
end
%==========================================================================

function plot_time_series_phase_cronologically(results,scrsz)
h_fig = figure(get_figure_h([scrsz(1)+round(scrsz(3)/4*2),scrsz(2)+52,...
    round(scrsz(3)/4),round(scrsz(4)-110)])); clf
setappdata(h_fig,'results',results)
setappdata(h_fig,'scrsz',scrsz)
set(h_fig,'Name',['Chrono modulation series' results.header.title])
uicontrol('Parent',h_fig, ...
	'BusyAction','cancel', ...
	'Callback',@CB_plot_prmSet, ...
	'Position',[65 1 60 15], ...
	'String','plot prmSet', ...
	'Style','pushbutton', ...
	'Tag','pushbutton_plot_prmSet', ...
	'UIContextMenu',[]);

result_names = fieldnames(results.data.presentation);
for q=1:length(result_names)
    eval(['series.' result_names{q} '= [results.data.presentation(:).' result_names{q} '];']);
end
plot_time_phase_series_spectral_lines(series)
%==========================================================================

function plot_time_series_phase_of_prmSet_averages(results,scrsz)

t = 'Function plot_time_series_phase_of_prmSet_averages not implemented yet'
