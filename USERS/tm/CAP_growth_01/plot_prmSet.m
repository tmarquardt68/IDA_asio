function plot_prmSet(results,prmSet,result,number,scrsz)
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
CAP_dur = length(result(number).CAP_unsuppr);
figure(get_figure_h([scrsz(1),scrsz(2)+round(scrsz(4)*1/4)-18,round(scrsz(3)/3.6),...
    round(scrsz(4)/4-40)]));clf
set(gcf,'Name',['CAP' results.header.title ' #' num2str(number)])

plot(linspace(0,1000/f_BT/8,CAP_dur),[result(number).CAP_unsuppr])
hold on
plot(linspace(0,2*1000/f_BT,2*8*CAP_dur),[result(number).modCAP; result(number).modCAP])

CAP_ampl2 = 2*(result(number).CAP_ampl-mean(result(number).CAP_ampl));
CAP_mod_curve = spline(5-80:10:75+80,[CAP_ampl2 CAP_ampl2 CAP_ampl2],linspace(0,80,8*CAP_dur));
plot([5:10:160]*2*1000/f_BT/160,[CAP_ampl2 CAP_ampl2],'kx')
plot(linspace(0,2*1000/f_BT,2*8*CAP_dur),[CAP_mod_curve CAP_mod_curve],'k-')

yLimit = get(gca,'YLim');
course = mean(yLimit) + diff(yLimit)/2 * result(number).course_BT_CAP;
plot(linspace(0,2*1000/f_BT,2*800), course,'k:'),

xlim([0 2/f_BT*1000]), grid on, zoom yon,
[~,idx]=max(course(1:end/2));
set(gca,'XTick',idx*100/3/1600+ [0 .25 .5 .75 1 ]/f_BT*1000,'FontSize',7)
set(gca,'XTickLabel',{'0' '.25' '.5' '.75' '1' '1.25' '1.5' '1.75' '2'})
title(['CAP modulation (prmSet:' num2str(prmSet) ') L_P=' num2str(l2) 'dB, L_B=' num2str(l_BT) 'dB    # ' num2str(number)])

% plot CM modulation patterns
figure(get_figure_h([scrsz(1),scrsz(2)+round(scrsz(4)*2/4)-18,round(scrsz(3)/3.6),...
    round(scrsz(4)/4-40)]));clf
set(gcf,'Name',['modCM patterns' results.header.title ' #' num2str(number)])
plot(x,result(number).modCM_2f1_f2_course,'b'),
hold on;
plot(x,20*log10(abs(result(number).l_unmodCM_2f1_f2)*ones(2/f_BT*fs,1)),'b--')
plot(x,result(number).modCM_f2_f1_course,'r'),
plot(x,20*log10(abs(result(number).l_unmodCM_f2_f1)*ones(2/f_BT*fs,1)),'r--')
% plot(x,(series(number).modCM_f2_course - 20*log10(abs(series(number).l_unmodCM_f2))),'k'), % shift to zero
% plot(x,zeros(2/f_BT*fs,1),'k--')
plot(x,result(number).modCM_2f2_f1_course,'y'),
plot(x,20*log10(abs(result(number).l_unmodCM_2f2_f1)*ones(2/f_BT*fs,1)),'y--')
plot(x,result(number).modCM_DPf1_course,'Color',[0 .9 0]),
plot(x,20*log10(abs(result(number).l_unmodCM_DPf1)*ones(2/f_BT*fs,1)),'g--','Color',[0 .9 0])
plot(x,result(number).modCM_DPf2_course,'w'),
plot(x,20*log10(abs(result(number).l_unmodCM_DPf2)*ones(2/f_BT*fs,1)),'w--')
plot(x,result(number).modCM_2f1_course,'m')
plot(x,20*log10(abs(result(number).l_unmodCM_2f1)*ones(2/f_BT*fs,1)),'m--')
plot(x,result(number).modCM_f2_course,'k'),
plot(x,20*log10(abs(result(number).l_unmodCM_f2)*ones(2/f_BT*fs,1)),'k--')

yLimit = get(gca,'YLim');
course = mean(yLimit) + diff(yLimit)/2 * result(number).course_BT;
plot(x, course,'k:'),
xlim([0 2/f_BT*1000]), grid on, zoom yon,
[~,idx]=max(course(1:end/2));
set(gca,'XTick',idx*100/3/1600+ [0 .25 .5 .75 1 ]/f_BT*1000,'FontSize',7)
set(gca,'XTickLabel',{'0' '.25' '.5' '.75' '1' '1.25' '1.5' '1.75' '2'})
set(gca,'Color',[.8 .8 .8])
title(['CM modulation (prmSet:' num2str(prmSet) ') L_P=' num2str(l2) 'dB, L_B=' num2str(l_BT) 'dB    # ' num2str(number)])

% plot OAE modulation patterns
figure(get_figure_h([scrsz(1),scrsz(2)+round(scrsz(4)*3/4)-18,...
    round(scrsz(3)/3.6),round(scrsz(4)/4-40)]));clf
set(gcf,'Name',['modOAE patterns' results.header.title ' #' num2str(number)])
plot(x,result(number).mod_2f1_f2_course,'b'),
hold on;
plot(x,20*log10(abs(result(number).l_unmod_2f1_f2)*ones(2/f_BT*fs,1)),'b--')
plot(x,result(number).mod_f2_f1_course,'r'),
plot(x,20*log10(abs(result(number).l_unmod_f2_f1)*ones(2/f_BT*fs,1)),'r--')
plot(x,result(number).mod_2f2_f1_course,'y'),
plot(x,20*log10(abs(result(number).l_unmod_2f2_f1)*ones(2/f_BT*fs,1)),'y--')
plot(x,result(number).modSF_DPf1_course,'Color',[0 .9 0]),
plot(x,20*log10(abs(result(number).l_unmod_SF_DPf1)*ones(2/f_BT*fs,1)),'g--','Color',[0 .9 0])
plot(x,result(number).modSF_DPf2_course,'w'),
plot(x,20*log10(abs(result(number).l_unmod_SF_DPf2)*ones(2/f_BT*fs,1)),'w--')
plot(x,result(number).modSF_f2_course,'k'),
plot(x,20*log10(abs(result(number).l_unmod_SF_f2)*ones(2/f_BT*fs,1)),'k--')
yLimit = get(gca,'YLim');
course = mean(yLimit) + diff(yLimit)/2 * result(number).course_BT;
plot(x, course,'k:'),
xlim([0 2/f_BT*1000]), grid on, zoom yon,
[~,idx]=max(course(1:end/2));
set(gca,'XTick',idx*100/3/1600+ [0 .25 .5 .75 1 ]/f_BT*1000,'FontSize',7)
set(gca,'XTickLabel',{'0' '.25' '.5' '.75' '1' '1.25' '1.5' '1.75' '2'})
set(gca,'Color',[.8 .8 .8])
title(['acoust. modulation (prmSet:' num2str(prmSet) ')  L_P=' num2str(l2) 'dB, L_B=' num2str(l_BT) 'dB    # ' num2str(number)])

% plot(bp_wave), hold on, plot(CAP_samples,zeros(length(CAP_samples),1),'.')

%% plot spectra
H_mod_DPodd = (result(number).H_modDP_1 - result(number).H_modDP_2)/2;
H_mod_DPeven = (result(number).H_modDP_1 + result(number).H_modDP_2)/2;
H_mod_DPall = (abs(result(number).H_modDP_1) + abs(result(number).H_modDP_2))/2;

% plot acoustic spectrum during modDP
figure(get_figure_h([scrsz(1)+round(scrsz(3)/3.6),scrsz(2)+round(scrsz(4)*3/4)-18,round(scrsz(3)/3.6),...
    round(scrsz(4)/4-40)]));clf
set(gcf,'Name',['OAE DPspectrum' results.header.title ' #' num2str(number)])
h=stem_60(linspace(0,fs/2-10,fs/20),20*log10(H_mod_DPall(1:fs/20,1)),'.');
set(h(1),'Color',[.5 .5 .5])
set(h(2),'Color',[.5 .5 .5])
hold on
stem_60((result(number).modlines_2f1_f2-1)*10,20*log10(abs(H_mod_DPodd(result(number).modlines_2f1_f2,1))),'b.')
stem_60((result(number).modlines_f2_f1-1)*10,20*log10(abs(H_mod_DPodd(result(number).modlines_f2_f1,1))),'r.')
stem_60((result(number).modlines_2f2_f1-1)*10,20*log10(abs(result(number).H_mod_2f2_f1)),'y.')
stem_60((result(number).modlines_f1-1)*10,20*log10(abs(result(number).H_mod_DPf1)),'g.')
stem_60((result(number).modlines_f2-1)*10,20*log10(abs(result(number).H_mod_DPf2)),'w.')
plot([f1 f2 f_BT],[l2+l1 l2 l_BT],'ko')
plot([f1-40 f2-40],[l2+l1+l_suppr l2+l_suppr],'ko')
plot([f1 f2],20*log10(abs([result(number).suppr_f1 result(number).suppr_f2])),'kx')
plot(f1,20*log10(abs(result(number).l_unmod_SF_DPf1)),'co')
plot(f2,20*log10(abs(result(number).l_unmod_SF_DPf2)),'co')
plot(f2-f1,20*log10(abs(result(number).l_unmod_f2_f1)),'co')
plot(2*f1-f2,20*log10(abs(result(number).l_unmod_2f1_f2)),'co')
plot(2*f2-f1,20*log10(abs(result(number).l_unmod_2f2_f1)),'co')
grid on,zoom on,axis([1 2.05*f2 -30 110]),set(gca,'FontSize',7),set(gca,'Color',[.8 .8 .8])
title(['Acoustic spectrum during modDP (prmSet:' num2str(prmSet) ')  L_P=' num2str(l2) 'dB, L_B=' num2str(l_BT) 'dB    # ' num2str(number)])

% plot acoustic spectrum during modSF
figure(get_figure_h([scrsz(1)+round(scrsz(3)/3.6*2),scrsz(2)+round(scrsz(4)*3/4)-18,round(scrsz(3)/3.6),...
    round(scrsz(4)/4-40)]));clf
set(gcf,'Name',['Acoustic SFspectum' results.header.title ' #' num2str(number)])
h=stem_60(linspace(0,fs/2-10,fs/20),20*log10(abs(result(number).H_mod_SF(1:fs/20,1))),'.');
set(h(1),'Color',[.5 .5 .5])
set(h(2),'Color',[.5 .5 .5])
hold on
stem_60((result(number).modlines_f2-1)*10,20*log10(abs(result(number).H_mod_f2)),'k.')
plot(f2,20*log10(abs(result(number).l_unmod_SF_f2)),'co')
grid on,zoom on,axis([1 2.05*f2 -30 110]),set(gca,'FontSize',7),set(gca,'Color',[.8 .8 .8])
title(['Acoustic spectrum during modSF (prmSet:' num2str(prmSet) ')  L_P=' num2str(l2) 'dB, L_B=' num2str(l_BT) 'dB    # ' num2str(number)])

% plot CM spectrum during modDP
figure(get_figure_h([scrsz(1)+round(scrsz(3)/3.6),scrsz(2)+round(scrsz(4)*2/4)-18,round(scrsz(3)/3.6),...
    round(scrsz(4)/4-40)]));clf
set(gcf,'Name',['CM DPspectum' results.header.title ' #' num2str(number)])
h = stem_60(linspace(0,fs/2-10,fs/20),20*log10(H_mod_DPall(1:fs/20,2)),'.');
set(h(1),'Color',[.5 .5 .5])
set(h(2),'Color',[.5 .5 .5])
hold on
stem_60((result(number).modlines_2f1_f2-1)*10,20*log10(abs(H_mod_DPodd(result(number).modlines_2f1_f2,2))),'b.')
stem_60((result(number).modlines_f2_f1-1)*10,20*log10(abs(H_mod_DPodd(result(number).modlines_f2_f1,2))),'r.')
stem_60((result(number).modlines_2f2_f1-1)*10,20*log10(abs(result(number).H_modCM_2f2_f1)),'y.')
stem_60((result(number).modlines_f1-1)*10,20*log10(abs(result(number).H_modCM_DPf1)),'g.')
stem_60((result(number).modlines_f2-1)*10,20*log10(abs(result(number).H_modCM_DPf2)),'w.')
stem_60((result(number).modlines_2f1-1)*10,20*log10(abs(H_mod_DPeven(result(number).modlines_2f1,2))),'m.')
plot(f2,20*log10(abs(result(number).l_unmodCM_DPf2)),'co')
plot(f1,20*log10(abs(result(number).l_unmodCM_DPf1)),'co')
plot(2*f1,20*log10(abs(result(number).l_unmodCM_2f1)),'co')
plot(f2-f1,20*log10(abs(result(number).l_unmodCM_f2_f1)),'co')
plot(2*f1-f2,20*log10(abs(result(number).l_unmodCM_2f1_f2)),'co')
plot(2*f2-f1,20*log10(abs(result(number).l_unmodCM_2f2_f1)),'co')
plot([f1 f2],20*log10(abs([result(number).supprCM_f1 result(number).supprCM_f2])),'kx')
grid on,zoom on,axis([1 2.05*f2 -30 60]),set(gca,'FontSize',7),set(gca,'Color',[.8 .8 .8])
title(['CM spectrum during modDP (prmSet:' num2str(prmSet) ')  L_P=' num2str(l2) 'dB, L_B=' num2str(l_BT) 'dB    # ' num2str(number)])

% plot CM spectrum during modSF
figure(get_figure_h([scrsz(1)+round(scrsz(3)/3.6*2),scrsz(2)+round(scrsz(4)*2/4)-18,...
    round(scrsz(3)/3.6),round(scrsz(4)/4-40)]));clf
set(gcf,'Name',['CM SFspectum' results.header.title ' #' num2str(number)])
h=stem_60(linspace(0,fs/2-10,fs/20),20*log10(abs(result(number).H_mod_SF(1:fs/20,2))),'.');
set(h(1),'Color',[.5 .5 .5])
set(h(2),'Color',[.5 .5 .5])
hold on
stem_60((result(number).modlines_f2-1)*10,20*log10(abs(result(number).H_mod_SF(result(number).modlines_f2,2))),'k.')
plot(f2,20*log10(abs(result(number).l_unmodCM_f2)),'co')
grid on,zoom on,axis([1 2.05*f2 -30 60]),set(gca,'FontSize',7),set(gca,'Color',[.8 .8 .8])
title(['CM spectrum during modSF (prmSet:' num2str(prmSet) ')  L_P=' num2str(l2) 'dB, L_B=' num2str(l_BT) 'dB    # ' num2str(number)])

