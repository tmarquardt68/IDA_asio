function plot_time_series_spectral_lines(series)

a=1;
% Plot modDP level
h_ax(a)=subplot('position',[0.0446    0.79    0.9429    0.18]);
plot(20*log10(abs(series.H_mod_2f1_f2(1,:))),'bo:'), hold on
plot(20*log10(abs(series.H_mod_2f1_f2(2,:))),'bx-')
plot(20*log10(abs(series.H_mod_2f1_f2(3,:))),'bx-','LineWidth',2),
plot(20*log10(abs(series.H_mod_2f1_f2(4,:))),'bx-')
plot(20*log10(abs(series.H_mod_2f1_f2(5,:))),'bo:')
plot(20*log10(abs(series.H_mod_f2_f1(1,:))),'ro:')
plot(20*log10(abs(series.H_mod_f2_f1(2,:))),'rx-')
plot(20*log10(abs(series.H_mod_f2_f1(3,:))),'ro-','LineWidth',2),
plot(20*log10(abs(series.H_mod_f2_f1(4,:))),'rx-')
plot(20*log10(abs(series.H_mod_f2_f1(5,:))),'ro:')
grid on,title modDP,set(gca,'FontSize',7); a=a+1;

% Plot modDP CM
h_ax(a)=subplot('position',[0.0446    0.60    0.9429    0.18]);
plot(20*log10(abs(series.H_modCM_2f1_f2(1,:))),'bo:'), hold on
plot(20*log10(abs(series.H_modCM_2f1_f2(2,:))),'bx-')
plot(20*log10(abs(series.H_modCM_2f1_f2(3,:))),'bx-','LineWidth',2),
plot(20*log10(abs(series.H_modCM_2f1_f2(4,:))),'bx-')
plot(20*log10(abs(series.H_modCM_2f1_f2(5,:))),'bo:')
plot(20*log10(abs(series.H_modCM_f2_f1(1,:))),'ro:')
plot(20*log10(abs(series.H_modCM_f2_f1(2,:))),'rx-')
plot(20*log10(abs(series.H_modCM_f2_f1(3,:))),'ro-','LineWidth',2),
plot(20*log10(abs(series.H_modCM_f2_f1(4,:))),'rx-')
plot(20*log10(abs(series.H_modCM_f2_f1(5,:))),'ro:')
grid on,title modDPCM,set(gca,'FontSize',7); a=a+1;

% Plot modSF level
h_ax(a)=subplot('position',[0.0446    0.41    0.9429    0.18]);
plot(20*log10(abs(series.H_mod_f2(1,:))),'ko:'), hold on
plot(20*log10(abs(series.H_mod_f2(2,:))),'kx-')
plot(20*log10(abs(series.H_mod_f2(3,:))),'k-','LineWidth',2),
plot(20*log10(abs(series.H_mod_f2(4,:))),'kx-')
plot(20*log10(abs(series.H_mod_f2(5,:))),'ko:')
plot(20*log10(abs(series.l_unmod_SF_f2)),'c.','LineWidth',2),
grid on,title modSF,set(gca,'FontSize',7); a=a+1;

% Plot modCM_f2 level
h_ax(a)=subplot('position',[0.0446   0.222   0.9429    0.18]);
plot(20*log10(abs(series.H_modCM_f2(1,:))),'ko:'), hold on
plot(20*log10(abs(series.H_modCM_f2(2,:))),'kx-')
plot(20*log10(abs(series.H_modCM_f2(3,:))),'kx-','LineWidth',2),
plot(20*log10(abs(series.H_modCM_f2(4,:))),'kx-')
plot(20*log10(abs(series.H_modCM_f2(5,:))),'ko:')
plot(20*log10(abs(series.l_unmodCM_f2)),'c.','LineWidth',2),
grid on,title modCMf2,set(gca,'FontSize',7); a=a+1;

% Plot CAPs
len_series = length(series.l_unmod_SF_f2);
h_ax(a)=subplot('position',[0.0446   0.03   0.9429    0.18]); 
%plot(reshape([series.CAP_ampl],8,len_series)'+repmat([0 1 2 3 4 5 6 7]* 0.05,len_series,1))
plot(reshape([series.CAP_ampl],8,len_series)')
grid on,title CAPs,set(gca,'FontSize',7)

xlabel('no of 6s-presentions'),zoom on
linkaxes(h_ax, 'x')
set(h_ax,'Color',[.8 .8 .8])
xlim(h_ax(1),[0 len_series])
zoom on
