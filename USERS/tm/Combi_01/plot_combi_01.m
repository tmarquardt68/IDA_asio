function plot_combi_01(filename,prmSet_to_plot)
% TO DO
% -

if ~exist('prmSet_to_plot','var'), prmSet_to_plot = 1;end

load(filename)

if isfield(results,'data')
    results = rmfield(results,'data');
end

monitor_settings = monitor_init_combi_01(results);

try
    result = monitor_combi_01(1,1,results,raw_1_1.wave,monitor_settings,1);
catch
    result = monitor_combi_01(1,1,results,raw_1.wave,monitor_settings,1);
end
no_prmSets = size(results.stimulus.original_parameter_table,1);
result_names = fieldnames(result);

for prmSet = no_prmSets:-1:1
    for q3=1:length(result_names)
        eval(['results.data.avg(' num2str(prmSet) ').' result_names{q3}...
            '=zeros(size(result.' result_names{q3} '));']);
    end

    % averaging
    for rpt=q-1:-1:1
        try
            eval(['wave = raw_' num2str(rpt) '_' num2str(prmSet) '.wave;']);
        catch
            eval(['wave = raw_' num2str(prmSet) '.wave;']);
        end

        result = monitor_combi_01(rpt,prmSet,results,wave,monitor_settings,0);
        results.data.sweeps(rpt,prmSet) = result;

        for q3=1:length(result_names)
            eval(['results.data.avg(prmSet).'  result_names{q3} '= results.data.avg(prmSet).' ...
                result_names{q3} '+ result.' result_names{q3} ';']);
        end
    end

    for q3=1:length(result_names)
        eval(['results.data.avg(prmSet).'  result_names{q3} '= results.data.avg(prmSet).' ...
            result_names{q3} '/(q-1);']);
    end
end
save(filename,'results','-append')

plot_prmSet(results,prmSet_to_plot,results.data.avg(prmSet_to_plot))
h_series = get_figure_h([1385 1 540 850]);
figure(h_series)
clf
set(h_series,'Name','Series')
plot(reshape([results.data.sweeps(:).CAP_ampl],8,q-1)'+repmat([0 1 2 3 4 5 6 7]* 0.05, q-1,1))

tmp =0;