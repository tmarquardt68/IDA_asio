function plot_combi_01(filename,prmSet_to_plot)

if ~exist('prmSet_to_plot','var'), prmSet_to_plot = 1;end

load(filename,'results')

if 1%~isfield(results,'data')
%     results = rmfield(results,'data');
    load(filename)
    monitor_settings = monitor_init_CAP_growth_01(results);

    result = monitor_CAP_growth_01(1,1,results,raw_1.wave,monitor_settings,0);
    results.data.avg(q2).mod_2f1_f2_course = zeros(length(result.mod_2f1_f2_course),1);
    results.data.avg(q2).mod_f2_f1_course = zeros(length(result.mod_f2_f1_course),1);
    results.data.avg(q2).modSF_f2 = zeros(length(result.modSF_f2),1);
    results.data.avg(q2).H_modCM_2f1_f2 = zeros(length(result.H_modCM_2f1_f2),1);
    results.data.avg(q2).H_mod_f2_f1 = zeros(length(result.H_mod_f2_f1),1);
    results.data.avg(q2).H_modCM_f2 = zeros(length(result.H_modCM_f2),1);
    results.data.avg(q2).CAP_unsuppr = zeros(length(result.CAP_unsuppr),1);
    results.data.avg(q2).modCAP = zeros(length(result.modCAP),1);
    results.data.course_BT(q2,:) = zeros(length(result.course_BT),1);

    for prmSet = q2:-1:1
        for rpt=q:-1:1
            eval(['wave = raw_' num2str(rpt) '_' num2str(prmSet) '.wave;']);
            result = monitor_CAP_growth_01(rpt,prmSet,results,wave,monitor_settings,0);

            % averaging
            results.data.avg(prmSet).mod_2f1_f2_course = results.data.avg(prmSet).mod_2f1_f2_course + result.mod_2f1_f2_course;
            results.data.avg(prmSet).mod_f2_f1_course = results.data.avg(prmSet).mod_f2_f1_course + result.mod_f2_f1_course;
            results.data.avg(prmSet).modSF_f2 = results.data.avg(prmSet).modSF_f2 + result.modSF_f2;
            results.data.avg(prmSet).H_modCM_2f1_f2 = results.data.avg(prmSet).H_modCM_2f1_f2 + result.H_modCM_2f1_f2;
            results.data.avg(prmSet).H_mod_f2_f1 = results.data.avg(prmSet).H_mod_f2_f1 + result.H_mod_f2_f1;
            results.data.avg(prmSet).H_modCM_f2 = results.data.avg(prmSet).H_modCM_f2 + result.H_modCM_f2;
            results.data.avg(prmSet).CAP_unsuppr = results.data.avg(prmSet).CAP_unsuppr + result.CAP_unsuppr;
            results.data.avg(prmSet).modCAP = results.data.avg(prmSet).modCAP + result.modCAP;

            results.data.sweeps(rpt,prmSet) = result;
        end
        results.data.course_BT(prmSet,:) = result.course_BT;
        results.data.avg(prmSet).mod_2f1_f2_course = results.data.avg(prmSet).mod_2f1_f2_course/q;
        results.data.avg(prmSet).mod_f2_f1_course = results.data.avg(prmSet).mod_f2_f1_course/q;
        results.data.avg(prmSet).modSF_f2 = results.data.avg(prmSet).modSF_f2/q;
        results.data.avg(prmSet).H_modCM_2f1_f2 = results.data.avg(prmSet).H_modCM_2f1_f2/q;
        results.data.avg(prmSet).H_mod_f2_f1 = results.data.avg(prmSet).H_mod_f2_f1/q;
        results.data.avg(prmSet).H_modCM_f2 = results.data.avg(prmSet).H_modCM_f2/q;
        results.data.avg(prmSet).CAP_unsuppr = results.data.avg(prmSet).CAP_unsuppr/q;
        results.data.avg(prmSet).modCAP = results.data.avg(prmSet).modCAP/q;
    end
    save(filename,'results','-append')
else 
    load(filename,'results')
end



%% plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_prmSet(results, results.data.avg(prmSet_to_plot))

% scrsz = get(0,'ScreenSize');
% fs = results.stimulus.sample_rate;
% f_BT = results.stimulus.parameters(2).values;
% f2 = results.stimulus.parameters(4).values;
% l2 = results.stimulus.parameters(6).values;
% l_BT = results.stimulus.parameters(8).values;
% x = linspace(0,2/f_BT*1000,2/f_BT*fs);
% 
% % plot acoustic modulation patterns
% position = [0,round(scrsz(4)*1/4)-18,...
%     round(scrsz(3)/3.6),round(scrsz(4)/4-40)];
% h_fig = findobj('Position',position);
% if (isempty(h_fig))
%     h_fig = figure;
%     set(h_fig,'Position',position)
% end
% figure(h_fig), clf,
% 
% hold off
% plot(x,results.data.avg(prmSet).mod_2f1_f2_course),
% hold on;
% plot(x,results.data.avg(prmSet).mod_f2_f1_course),
% plot(x,results.data.avg(prmSet).modSF_f2),
% 
% yLimit = get(gca,'YLim');
% course = mean(yLimit) + diff(yLimit)/2 * results.data.course_BT(prmSet,:);
% plot(x, course,'k:'),
% xlim([0 2/f_BT*1000]), grid on, zoom on,
% legend({'2f2-f1' 'f2-f1' 'f2'})
% 
% % plot CM modulation patterns
% position = [round(scrsz(3)/3.6),round(scrsz(4)*1/4)-18,...
%     round(scrsz(3)/3.6),round(scrsz(4)/4-40)];
% h_fig = findobj('Position',position);
% if (isempty(h_fig))
%     h_fig = figure;
%     set(h_fig,'Position',position)
% end
% figure(h_fig), clf,
% 
% hold off
% plot(x,results.data.avg(prmSet).H_modCM_2f1_f2),
% hold on;
% plot(x,results.data.avg(prmSet).H_mod_f2_f1),
% 
% yLimit = get(gca,'YLim');
% course = mean(yLimit) + 0.8*diff(yLimit)/2 * ...
%     (results.data.avg(prmSet).H_modCM_f2-mean(results.data.avg(prmSet).H_modCM_f2))/...
%     max((results.data.avg(prmSet).H_modCM_f2-mean(results.data.avg(prmSet).H_modCM_f2)));
% plot(x,course),
% course = mean(yLimit) + diff(yLimit)/2 * results.data.course_BT(prmSet,:);
% plot(x, course,'k:'),
% xlim([0 2/f_BT*1000]), grid on, zoom on,
% legend({'2f2-f1' 'f2-f1' 'f2 (a.u.)'})
% 
% 
% % plot CAPs
% CAP_dur = length(results.data.avg(prmSet).CAP_unsuppr);
% position = [round(scrsz(3)/3.6*2),round(scrsz(4)*3/4)-18,...
%     round(scrsz(3)/3.6),round(scrsz(4)/4-40)];
% h_fig = findobj('Position',position);
% if (isempty(h_fig))
%     h_fig = figure;
%     set(h_fig,'Position',position)
% end
% figure(h_fig), clf,
% plot(linspace(0,10,CAP_dur),results.data.avg(prmSet).CAP_unsuppr)
% hold on
% plot(linspace(10,90,8*CAP_dur),results.data.avg(prmSet).modCAP)
% grid on, % axis([1 24000 -30 70])
% title(['f_BT=' num2str(f_BT) 'Hz @ ' num2str(l_BT) ' dB;  f2=' num2str(f2) 'HZ @ ' num2str(l2) ' dB.'])
% 
% % plot(bp_wave), hold on, plot(CAP_samples,zeros(length(CAP_samples),1),'.')
% 
