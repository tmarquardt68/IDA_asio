function results = plot_combi_01(results,filename)
% TO DO
% -

monitor_no = 3; % choose monitor


scrsz = get(0,'MonitorPositions');
if monitor_no > size(scrsz,1)
    monitor_no = 1; 
end
data_path = [fileparts(which(mfilename)) '\..\..\..\DATA\'];

if isempty(results) && exist('filename','Var') % Data field is only filled the first time
    load([data_path filename],'results')
else
    error('give either $resutls or $filename as parameter!')
end
if ~isfield(results,'data') % plot combi runs the first time (monitor_combi does only save raw recordings)
    no_prmSets = size(results.stimulus.original_parameter_table,1);
    monitor_settings = monitor_init_combi_01(results);
    load([data_path results.header.title])
    results = rmfield(results,'data');
    interval_order = results.stimulus.interval_order';
    interval_order = interval_order(1:end);
    for q = 1:length(interval_order)
        try
            eval(['wave = raw_' num2str(q) '.wave;']);
        catch % this catch is just for spme older filename convention
            eval(['wave = raw_' num2str(ceil(q/no_prmSets)) '_' num2str(mod(q-1,no_prmSets)+1) '.wave;']);
        end
        results.data.presentation(q) = monitor_combi_01([],interval_order(q),results,wave,monitor_settings,0);
    end
    save([data_path filename],'results','-append')
end

plot_time_series_spectral_line_cronologically(results,scrsz(monitor_no,:));
setappdata(gcf,'initial_call',1)
CB_plot_prmSet

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot time series of relevant spectral line levels

function plot_time_series_spectral_line_cronologically(results,scrsz)
series.H_mod_2f1_f2 = [results.data.presentation(:).H_mod_2f1_f2];
series.H_mod_f2_f1 = [results.data.presentation(:).H_mod_f2_f1];
series.H_modCM_f2_f1 = [results.data.presentation(:).H_modCM_f2_f1];
series.H_modCM_2f1_f2 = [results.data.presentation(:).H_modCM_2f1_f2];
series.H_mod_f2 = [results.data.presentation(:).H_mod_f2];
series.l_SF = [results.data.presentation(:).l_SF];
series.H_modCM_f2 = [results.data.presentation(:).H_modCM_f2];
series.l_unmodCM_f2 = [results.data.presentation(:).l_unmodCM_f2];
series.CAP_ampl = [results.data.presentation(:).CAP_ampl];

h_fig = figure(get_figure_h([scrsz(1)+round(scrsz(3)/4*3),scrsz(2)+40,...
    round(scrsz(3)/4),round(scrsz(4)-100)])); clf
setappdata(h_fig,'results',results)
setappdata(h_fig,'series',series)
setappdata(h_fig,'scrsz',scrsz)
set(h_fig,'Name',['Chrono series' results.header.title])
uicontrol('Parent',h_fig, ...
	'BusyAction','cancel', ...
	'Callback',@CB_toggle_series, ...
	'Position',[1 1 60 15], ...
	'String','plot avg', ...
	'Style','pushbutton', ...
	'Tag','pushbutton_toggle', ...
	'UIContextMenu',[]);
uicontrol('Parent',h_fig, ...
	'BusyAction','cancel', ...
	'Callback',@CB_plot_prmSet, ...
	'Position',[65 1 60 15], ...
	'String','plot prmSet', ...
	'Style','pushbutton', ...
	'Tag','pushbutton_plot_prmSet', ...
	'UIContextMenu',[]);
plot_time_series_spectral_lines(series)
%==========================================================================

function plot_time_series_spectral_lines_of_prmSet_averages(results,scrsz)
no_prmSets = size(results.stimulus.original_parameter_table,1);
result_names = fieldnames(results.data.presentation);
interval_order = results.stimulus.interval_order';
interval_order = interval_order(1:end);

for prmSet = no_prmSets:-1:1
    for q=1:length(result_names)
        eval(['avg(' num2str(prmSet) ').' result_names{q}...
            '=zeros(size(results.data.presentation(1).' result_names{q} '));']);
    end
    n=0;
    for q2=1:length(interval_order)
        if interval_order(q2) == prmSet
            n=n+1;
            for q=1:length(result_names)
                eval(['avg(prmSet).'  result_names{q} '= avg(prmSet).' ...
                    result_names{q} '+ (results.data.presentation(q2).' result_names{q} ');']);
            end
        end
    end


    for q=1:length(result_names)
        eval(['avg(prmSet).'  result_names{q} '= avg(prmSet).' ...
            result_names{q} '/n;']);
    end
end

series.H_mod_2f1_f2 = [avg(:).H_mod_2f1_f2];
series.H_mod_f2_f1 = [avg(:).H_mod_f2_f1];
series.H_modCM_f2_f1 = [avg(:).H_modCM_f2_f1];
series.H_modCM_2f1_f2 = [avg(:).H_modCM_2f1_f2];
series.H_mod_f2 = [avg(:).H_mod_f2];
series.l_SF = [avg(:).l_SF];
series.H_modCM_f2 = [avg(:).H_modCM_f2];
series.l_unmodCM_f2 = [avg(:).l_unmodCM_f2];
series.CAP_ampl = [avg(:).CAP_ampl];

h_fig = figure(get_figure_h([scrsz(1)+round(scrsz(3)/4*3),scrsz(2)+40,...
    round(scrsz(3)/4),round(scrsz(4)-100)])); clf
results.data.avg = avg;
setappdata(h_fig,'results',results)
setappdata(h_fig,'series',series)
setappdata(h_fig,'scrsz',scrsz)
set(h_fig,'Name',['Avg series' results.header.title])
uicontrol('Parent',h_fig, ...
	'BusyAction','cancel', ...
	'Callback',@CB_toggle_series, ...
	'Position',[1 1 60 10], ...
	'String','plot chrono', ...
	'Style','pushbutton', ...
	'Tag','pushbutton_toggle', ...
	'UIContextMenu',[]);
uicontrol('Parent',h_fig, ...
	'BusyAction','cancel', ...
	'Callback',@CB_plot_prmSet, ...
	'Position',[65 1 60 10], ...
	'String','plot prmSet', ...
	'Style','pushbutton', ...
	'Tag','pushbutton_plot_prmSet', ...
	'UIContextMenu',[]);
plot_time_series_spectral_lines(series)
%==========================================================================

function CB_plot_prmSet(callingObj,~)
% get number from mouse click position
if getappdata(gcf,'initial_call')
    number = [];
    setappdata(gcf,'initial_call',0)
else
    [x, ~] = ginput(1);
    number = round(x);
    zoom on
end

results = getappdata(gcf,'results');
scrsz = getappdata(gcf,'scrsz');
kindOfPlot = get(gcf,'Name');

if strcmp(kindOfPlot(1:3),'Avg')
    if isempty(number)
        number = round(size(results.stimulus.original_parameter_table,1)/1);
    end
    if number>=1 || number<=size(results.stimulus.original_parameter_table,1)
        plot_prmSet(results,number,results.data.avg,number,scrsz)
    end
else
    interval_order = results.stimulus.interval_order';
    if isempty(number)
        number = round(length(interval_order)/2);
    end
    prmSet = interval_order(number);
    if number>=1 || number<=length(interval_order)
        plot_prmSet(results,prmSet,results.data.presentation,number,scrsz)
    end
end
%==========================================================================

function CB_toggle_series(callingObj,~)
results = getappdata(gcbf,'results');
scrsz = getappdata(gcbf,'scrsz');
kindOfPlot = get(gcbf,'Name');

if strcmp(kindOfPlot(1:3),'Avg')
    plot_time_series_spectral_line_cronologically(results,scrsz)
    set(findobj('Tag','pushbutton_toggle'),'String','plot_avg')
else
    plot_time_series_spectral_lines_of_prmSet_averages(results,scrsz);
    set(findobj('Tag','pushbutton_toggle'),'String','plot_chrono')
end
