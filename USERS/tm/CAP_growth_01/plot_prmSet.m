function plot_prmSet(results,prmSet,series,number,scrsz)
% TO DO:
% - plot BT waveform in CAP modualtion figure, scale spline from buttom of
% axis = 0
% - cyclic fit of spline
% - plot also CM modulation time course of 1st harmonicof f2 abd it's
% intermoduilations procusts with f2. 
% - plot also of modualtion time course of the f1's and f2's CM and SFOAE.

fs = results.stimulus.sample_rate;
f1 = results.stimulus.original_parameter_table(prmSet,3);
f2 = results.stimulus.original_parameter_table(prmSet,4);
f_BT = results.stimulus.original_parameter_table(prmSet,2);
l1 = results.stimulus.original_parameter_table(prmSet,5);
l2 = results.stimulus.original_parameter_table(prmSet,6);
l_suppr = results.stimulus.original_parameter_table(prmSet,7);
l_BT = results.stimulus.original_parameter_table(prmSet,8);
x = linspace(0,2*1000/f_BT,2/f_BT*fs);

%% plot time courses

% plot CAPs
CAP_dur = length(series(number).CAP_unsuppr);
figure(get_figure_h([scrsz(1),scrsz(2)+round(scrsz(4)*1/4)-18,round(scrsz(3)/3.6),...
    round(scrsz(4)/4-40)]));clf
set(gcf,'Name',['CAP' results.header.title ' #' num2str(number)])

plot(linspace(0,1000/f_BT/8,CAP_dur),[series(number).CAP_unsuppr])
hold on
plot(linspace(0,2*1000/f_BT,2*8*CAP_dur),[series(number).modCAP; series(number).modCAP])

CAP_ampl2 = 2*(series(number).CAP_ampl-mean(series(number).CAP_ampl));
CAP_mod_curve = spline(5-80:10:75+80,[CAP_ampl2 CAP_ampl2 CAP_ampl2],linspace(0,80,8*CAP_dur));
plot([5:10:160]*2*1000/f_BT/160,[CAP_ampl2 CAP_ampl2],'kx')
plot(linspace(0,2*1000/f_BT,2*8*CAP_dur),[CAP_mod_curve CAP_mod_curve],'k-')

yLimit = get(gca,'YLim');
course = mean(yLimit) + diff(yLimit)/2 * series(number).course_BT_CAP;
plot(linspace(0,2*1000/f_BT,2*800), course,'k:'),

xlim([0 2/f_BT*1000]), grid on, zoom on,
[~,idx]=max(course(1:end/2));
set(gca,'XTick',idx*100/3/1600+ [0 .25 .5 .75 1 ]/f_BT*1000,'FontSize',7)
set(gca,'XTickLabel',{'0' '.25' '.5' '.75' '1' '1.25' '1.5' '1.75' '2'})
title(['CAP: L_P=' num2str(l2) 'dB, L_B=' num2str(l_BT) 'dB    # ' num2str(number)])

% plot CM modulation patterns
figure(get_figure_h([scrsz(1),scrsz(2)+round(scrsz(4)*2/4)-18,round(scrsz(3)/3.6),...
    round(scrsz(4)/4-40)]));clf
set(gcf,'Name',['modCM patterns' results.header.title ' #' num2str(number)])
plot(x,series(number).modCM_2f1_f2_course,'b'),
hold on;
plot(x,20*log10(abs(series(number).l_unmodCM_2f1_f2)*ones(2/f_BT*fs,1)),'b--')
plot(x,series(number).modCM_f2_f1_course,'r'),
plot(x,20*log10(abs(series(number).l_unmodCM_f2_f1)*ones(2/f_BT*fs,1)),'r--')
% plot(x,(series(number).modCM_f2_course - 20*log10(abs(series(number).l_unmodCM_f2))),'k'), % shift to zero
% plot(x,zeros(2/f_BT*fs,1),'k--')
plot(x,(series(number).modCM_f2_course),'k'),
plot(x,20*log10(abs(series(number).l_unmodCM_f2)*ones(2/f_BT*fs,1)),'k--')

yLimit = get(gca,'YLim');
course = mean(yLimit) + diff(yLimit)/2 * series(number).course_BT;
plot(x, course,'k:'),
xlim([0 2/f_BT*1000]), grid on, zoom on,
[~,idx]=max(course(1:end/2));
set(gca,'XTick',idx*100/3/1600+ [0 .25 .5 .75 1 ]/f_BT*1000,'FontSize',7)
set(gca,'XTickLabel',{'0' '.25' '.5' '.75' '1' '1.25' '1.5' '1.75' '2'})
title(['CM modulation: L_P=' num2str(l2) 'dB, L_B=' num2str(l_BT) 'dB    # ' num2str(number)])

% plot OAE modulation patterns
figure(get_figure_h([scrsz(1),scrsz(2)+round(scrsz(4)*3/4)-18,...
    round(scrsz(3)/3.6),round(scrsz(4)/4-40)]));clf
set(gcf,'Name',['modOAE patterns' results.header.title ' #' num2str(number)])
plot(x,series(number).mod_2f1_f2_course,'b'),
hold on;
plot(x,20*log10(abs(series(number).l_2f1_f2)*ones(2/f_BT*fs,1)),'b--')
plot(x,series(number).mod_f2_f1_course,'r'),
plot(x,20*log10(abs(series(number).l_f2_f1)*ones(2/f_BT*fs,1)),'r--')
plot(x,series(number).modSF_f2,'k'),
plot(x,20*log10(abs(series(number).l_SF)*ones(2/f_BT*fs,1)),'k--')
yLimit = get(gca,'YLim');
course = mean(yLimit) + diff(yLimit)/2 * series(number).course_BT;
plot(x, course,'k:'),
xlim([0 2/f_BT*1000]), grid on, zoom on,
[~,idx]=max(course(1:end/2));
set(gca,'XTick',idx*100/3/1600+ [0 .25 .5 .75 1 ]/f_BT*1000,'FontSize',7)
set(gca,'XTickLabel',{'0' '.25' '.5' '.75' '1' '1.25' '1.5' '1.75' '2'})
title(['acoust. modulation: L_P=' num2str(l2) 'dB, L_B=' num2str(l_BT) 'dB    # ' num2str(number)])

% plot(bp_wave), hold on, plot(CAP_samples,zeros(length(CAP_samples),1),'.')

%% plot spectra
H_mod_DPodd = (series(number).H_modDP_1 - series(number).H_modDP_2)/2;
H_mod_DPall = (abs(series(number).H_modDP_1) + abs(series(number).H_modDP_2))/2;

% plot acoustic spectrum during modDP
figure(get_figure_h([scrsz(1)+round(scrsz(3)/3.6),scrsz(2)+round(scrsz(4)*3/4)-18,round(scrsz(3)/3.6),...
    round(scrsz(4)/4-40)]));clf
set(gcf,'Name',['OAE DPspectrum' results.header.title ' #' num2str(number)])
h=stem_60(linspace(0,fs/2-10,fs/20),20*log10(H_mod_DPall(1:fs/20,1)),'.');
set(h(1),'Color',[.5 .5 .5])
set(h(2),'Color',[.5 .5 .5])
hold on
stem_60((series(number).modlines_2f1_f2-1)*10,20*log10(abs(H_mod_DPodd(series(number).modlines_2f1_f2,1))),'b.')
stem_60((series(number).modlines_f2_f1-1)*10,20*log10(abs(H_mod_DPodd(series(number).modlines_f2_f1,1))),'r.')
plot([f1 f2 f_BT],[l2+l1 l2 l_BT],'ko')
plot([f1-40 f2-40],[l2+l1+l_suppr l2+l_suppr],'ko')
plot([f1-40 f2-40],[20*log10(abs(series(number).l1_suppr)) ...
    20*log10(abs(series(number).l2_suppr))],'kx')
plot(f2,20*log10(abs(series(number).l_SF_DPf2)),'mo')
plot(f1,20*log10(abs(series(number).l_SF_DPf1)),'mo')
plot(f2,20*log10(abs(series(number).l_SF_modDPf2)),'go')
plot(f1,20*log10(abs(series(number).l_SF_modDPf1)),'go')
grid on,zoom on,title(['Acoustic spectrum during modDP: L_P=' num2str(l2) 'dB, L_B=' num2str(l_BT) 'dB    # ' num2str(number)])
axis([1 2.05*f2 -30 110]),set(gca,'FontSize',7)

% plot acoustic spectrum during modSF
figure(get_figure_h([scrsz(1)+round(scrsz(3)/3.6*2),scrsz(2)+round(scrsz(4)*3/4)-18,round(scrsz(3)/3.6),...
    round(scrsz(4)/4-40)]));clf
set(gcf,'Name',['Acoustic SFspectum' results.header.title ' #' num2str(number)])
h=stem_60(linspace(0,fs/2-10,fs/20),20*log10(abs(series(number).H_mod_SF(1:fs/20,1))),'.');
set(h(1),'Color',[.5 .5 .5])
set(h(2),'Color',[.5 .5 .5])
hold on
stem_60((series(number).modlines_f2([1 2 4 5])-1)*10,20*log10(abs(series(number).H_mod_SF(series(number).modlines_f2([1 2 4 5]) ,1))),'k.')
plot(f2,20*log10(abs(series(number).l_SF)),'mo')
plot(f2,20*log10(abs(series(number).l_modSF)),'ko')
grid on,zoom on,title(['Acoustic spectrum during modSF: L_P=' num2str(l2) 'dB, L_B=' num2str(l_BT) 'dB    # ' num2str(number)])
axis([1 2.05*f2 -30 110]),set(gca,'FontSize',7)

% plot CM spectrum during modDP
figure(get_figure_h([scrsz(1)+round(scrsz(3)/3.6),scrsz(2)+round(scrsz(4)*2/4)-18,round(scrsz(3)/3.6),...
    round(scrsz(4)/4-40)]));clf
set(gcf,'Name',['CM DPspectum' results.header.title ' #' num2str(number)])
h = stem_60(linspace(0,fs/2-10,fs/20),20*log10(H_mod_DPall(1:fs/20,2)),'.');
set(h(1),'Color',[.5 .5 .5])
set(h(2),'Color',[.5 .5 .5])
hold on
stem_60((series(number).modlines_2f1_f2-1)*10,20*log10(abs(H_mod_DPodd(series(number).modlines_2f1_f2,2))),'b.')
stem_60((series(number).modlines_f2_f1-1)*10,20*log10(abs(H_mod_DPodd(series(number).modlines_f2_f1,2))),'r.')
plot(f2,20*log10(abs(series(number).l_unmodCM_DPf2)),'mo')
plot(f1,20*log10(abs(series(number).l_unmodCM_DPf1)),'mo')
plot(f1,20*log10(abs(series(number).supprCM_l1)),'kx')
grid on,zoom on,title(['CM spectrum during modDP: L_P=' num2str(l2) 'dB, L_B=' num2str(l_BT) 'dB    # ' num2str(number)])
axis([1 2.05*f2 -30 60]),set(gca,'FontSize',7)

% plot CM spectrum during modSF
figure(get_figure_h([scrsz(1)+round(scrsz(3)/3.6*2),scrsz(2)+round(scrsz(4)*2/4)-18,...
    round(scrsz(3)/3.6),round(scrsz(4)/4-40)]));clf
set(gcf,'Name',['CM SFspectum' results.header.title ' #' num2str(number)])
h=stem_60(linspace(0,fs/2-10,fs/20),20*log10(abs(series(number).H_mod_SF(1:fs/20,2))),'.');
set(h(1),'Color',[.5 .5 .5])
set(h(2),'Color',[.5 .5 .5])
hold on
stem_60((series(number).modlines_f2-1)*10,20*log10(abs(series(number).H_mod_SF(series(number).modlines_f2,2))),'k.')
plot(f2,20*log10(abs(series(number).l_unmodCM_f2)),'mo')
plot(f2,20*log10(abs(series(number).supprCM_l2)),'kx')
grid on,zoom on,title(['CM spectrum during modSF: L_P=' num2str(l2) 'dB, L_B=' num2str(l_BT) 'dB    # ' num2str(number)])
axis([1 2.05*f2 -30 60]),set(gca,'FontSize',7)

