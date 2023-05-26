function plot_CAP_01(results)

[h_subplots, ~] = gensubpl('hblcfg',5, 12,1, 1);
static_parameters = results.stimulus.stimuli_parameters.stimuli.specs.param;
parameter_table = results.stimulus.original_parameter_table(results.stimulus.interval_order,:);
sample_rate = results.stimulus.stimuli_parameters.sample_rate;

for q2=1:length(results.raw)
    eCoch_wave=results.raw(q2).wave(:,2);

    %%% plot results (quick fix)
    curr_parameters = parameter_table(q2,2:size(parameter_table,2));
    plot_length = curr_parameters(2)+20 /1000*sample_rate; % samples
    %period = 2*static_parameters(5).value+curr_parameters(2)+curr_parameters(4);
    n_reps = static_parameters(6).value;
    period_len = size(eCoch_wave,1)/n_reps;
    n_reps = n_reps-2;
    eCoch_wave_periods = reshape(eCoch_wave(2001:n_reps*period_len+2000),period_len,[]);
    subplot(h_subplots(q2))
    hold off
    plot(linspace(0,plot_length*1000/sample_rate,plot_length), ... CM
        (sum(eCoch_wave_periods(1:plot_length,1:2:end),2) - ...
        sum(eCoch_wave_periods(1:plot_length,2:2:end),2))/n_reps,'b')
    hold on
    neg_peak(q2) = min(sum(eCoch_wave_periods(1:plot_length,:),2)/n_reps);
    pos_peak(q2) = max(sum(eCoch_wave_periods(1:plot_length,:),2)/n_reps);
    plot(linspace(0,plot_length*1000/sample_rate,plot_length),... CAP
        sum(eCoch_wave_periods(1:plot_length,:),2)/n_reps,'r-')
    plot(linspace(0,plot_length*1000/sample_rate,plot_length),...
        neg_peak(q2)*ones(1,size(eCoch_wave_periods(1:plot_length,:),1)),'r--')
    plot(linspace(0,plot_length*1000/sample_rate,plot_length),...
        pos_peak(q2)*ones(1,size(eCoch_wave_periods(1:plot_length,:),1)),'r--')
    zoom on, grid on,
end

for q=1:length(h_subplots)
    subplot(h_subplots(q))
    ylim([min(neg_peak)*1.1, max(pos_peak)*1.1])
end