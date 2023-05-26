function plot_prmSet(results, prmSet, result)
% TO DO:
% - legnds are wrong!
% - plot average spectrum
% - plot BT waveform in CAP modualtion figure, scale spline from buttom of
% axis = 0
% - cyclic fit of spline

scrsz = get(0,'ScreenSize');
fs = results.stimulus.sample_rate;
f_BT = results.stimulus.original_parameter_table(prmSet,2);
l2 = results.stimulus.original_parameter_table(prmSet,6);
l_BT = results.stimulus.original_parameter_table(prmSet,8);

% plot CAPs
CAP_dur = length(result.CAP_unsuppr);
figure(get_figure_h([0,round(scrsz(4)*1/4)-18,round(scrsz(3)/3.6),...
    round(scrsz(4)/4-40)]));clf
plot(linspace(0,10,CAP_dur),result.CAP_unsuppr)
hold on
plot(linspace(10,90,8*CAP_dur),result.modCAP)
for q =1:8
    CAP_ampl(q) = max(result.modCAP((q-1)*CAP_dur+1:q*CAP_dur)) - ...
        min(result.modCAP((q-1)*CAP_dur+1:q*CAP_dur));
end
CAP_ampl2 = 2*(CAP_ampl-mean(CAP_ampl));
CAP_mod_curve = spline(5:10:95,[CAP_ampl2(end) CAP_ampl2 CAP_ampl2(1)],linspace(10,90,8*CAP_dur));
plot(15:10:85,2*(CAP_ampl-mean(CAP_ampl)),'kx')
plot(linspace(10,90,8*CAP_dur),CAP_mod_curve,'k-')
grid on, title(['CAP - L_P=' num2str(l2) 'dB, L_B=' num2str(l_BT) 'dB']) % axis([1 24000 -30 70])

%% plot time courses
x = linspace(0,2/f_BT*1000,2/f_BT*fs);

% plot CM modulation patterns
figure(get_figure_h([0,round(scrsz(4)*2/4)-18,round(scrsz(3)/3.6),...
    round(scrsz(4)/4-40)]));clf
plot(x,result.modCM_2f1_f2_course,'b'),
hold on;
plot(x,20*log10(abs(result.l_unmodCM_2f1_f2)*ones(2/f_BT*fs,1)),'b--')
plot(x,result.modCM_f2_f1_course,'r'),
plot(x,20*log10(abs(result.l_unmodCM_f2_f1)*ones(2/f_BT*fs,1)),'r--')
plot(x,10*(result.modCM_f2_course - 20*log10(abs(result.l_unmodCM_f2))),'k'),
plot(x,zeros(2/f_BT*fs,1),'k--')
yLimit = get(gca,'YLim');
course = mean(yLimit) + diff(yLimit)/2 * result.course_BT;
plot(x, course,'k:'),
xlim([0 2/f_BT*1000]), grid on, zoom on,
legend({'2f2-f1' 'f2-f1' '10xf2'}), title CM

% plot acoustic modulation patterns
figure(get_figure_h([0,round(scrsz(4)*3/4)-18,...
    round(scrsz(3)/3.6),round(scrsz(4)/4-40)]));clf
plot(x,result.mod_2f1_f2_course,'b'),
hold on;
plot(x,20*log10(abs(result.l_2f1_f2)*ones(2/f_BT*fs,1)),'b--')
plot(x,result.mod_f2_f1_course,'r'),
plot(x,20*log10(abs(result.l_f2_f1)*ones(2/f_BT*fs,1)),'r--')
plot(x,result.modSF_f2,'k'),
plot(x,20*log10(abs(result.l_SF)*ones(2/f_BT*fs,1)),'k--')

l = length(results.H_mic);
H = zeros(l,1);
H(result.modlines_f2_f1) = result.H_mod_f2_f1;
course = abs(ifft(H));
plot(x,20*log10(l * course(1:2/f_BT*fs)),'m');

yLimit = get(gca,'YLim');
course = mean(yLimit) + diff(yLimit)/2 * result.course_BT;
plot(x, course,'k:'),
xlim([0 2/f_BT*1000]), grid on, zoom on,
legend({'2f2-f1' 'f2-f1' 'f2'}), title OAE

% plot(bp_wave), hold on, plot(CAP_samples,zeros(length(CAP_samples),1),'.')
