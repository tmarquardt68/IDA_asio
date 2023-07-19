function plot_time_magn_series_spectral_lines(series)

a=1;
% Plot modDP level
h_ax(a)=subplot('position',[0.0446    0.79    0.9429    0.18]);
plot(20*log10(abs(series.H_mod_2f2_f1(1,:))),'y.:'), hold on
plot(20*log10(abs(series.H_mod_2f2_f1(2,:))),'yx-')
plot(20*log10(abs(series.H_mod_2f2_f1(3,:))),'y-','LineWidth',3),
plot(20*log10(abs(series.H_mod_2f2_f1(4,:))),'yx-')
plot(20*log10(abs(series.H_mod_2f2_f1(5,:))),'y.:')
plot(20*log10(abs(series.l_unmod_2f2_f1)),'c.','MarkerEdgeColor',[.9 .9 .5])
plot(20*log10(abs(series.H_mod_2f1(1,:))),'m.:'), hold on
plot(20*log10(abs(series.H_mod_2f1(2,:))),'mx-')
plot(20*log10(abs(series.H_mod_2f1(3,:))),'m-','LineWidth',3),
plot(20*log10(abs(series.H_mod_2f1(4,:))),'mx-')
plot(20*log10(abs(series.H_mod_2f1(5,:))),'m.:') 
% plot(20*log10(abs(series.l_unmod_2f1)),'c.','MarkerEdgeColor',[1 .6 1]) %curently not a filed in $results
plot(20*log10(abs(series.H_mod_2f1_f2(1,:))),'b.:'), hold on
plot(20*log10(abs(series.H_mod_2f1_f2(2,:))),'bx-')
plot(20*log10(abs(series.H_mod_2f1_f2(3,:))),'b-','LineWidth',3),
plot(20*log10(abs(series.H_mod_2f1_f2(4,:))),'bx-')
plot(20*log10(abs(series.H_mod_2f1_f2(5,:))),'b.:')
plot(20*log10(abs(series.l_unmod_2f1_f2)),'c.','MarkerEdgeColor',[.6 .6 1])
plot(20*log10(abs(series.H_mod_f2_f1(1,:))),'r.:')
plot(20*log10(abs(series.H_mod_f2_f1(2,:))),'rx-')
plot(20*log10(abs(series.H_mod_f2_f1(3,:))),'r-','LineWidth',3),
plot(20*log10(abs(series.H_mod_f2_f1(4,:))),'rx-')
plot(20*log10(abs(series.H_mod_f2_f1(5,:))),'r.:')
plot(20*log10(abs(series.l_unmod_f2_f1)),'c.','MarkerEdgeColor',[1 .6 .6])
grid on,title modDP,set(gca,'FontSize',7); a=a+1;

% Plot modDP CM
h_ax(a)=subplot('position',[0.0446    0.60    0.9429    0.18]);
plot(20*log10(abs(series.H_modCM_2f2_f1(1,:))),'y.:'), hold on
plot(20*log10(abs(series.H_modCM_2f2_f1(2,:))),'yx-')
plot(20*log10(abs(series.H_modCM_2f2_f1(3,:))),'y-','LineWidth',3),
plot(20*log10(abs(series.H_modCM_2f2_f1(4,:))),'yx-')
plot(20*log10(abs(series.H_modCM_2f2_f1(5,:))),'y.:')
plot(20*log10(abs(series.l_unmodCM_2f2_f1)),'c.','MarkerEdgeColor',[.9 .9 .5])
plot(20*log10(abs(series.H_modCM_2f1(1,:))),'m.:'), hold on
plot(20*log10(abs(series.H_modCM_2f1(2,:))),'mx-')
plot(20*log10(abs(series.H_modCM_2f1(3,:))),'m-','LineWidth',3),
plot(20*log10(abs(series.H_modCM_2f1(4,:))),'mx-')
plot(20*log10(abs(series.H_modCM_2f1(5,:))),'m.:')
plot(20*log10(abs(series.l_unmodCM_2f1)),'c.','MarkerEdgeColor',[1 .6 1])
plot(20*log10(abs(series.H_modCM_2f1_f2(1,:))),'b.:'), hold on
plot(20*log10(abs(series.H_modCM_2f1_f2(2,:))),'bx-')
plot(20*log10(abs(series.H_modCM_2f1_f2(3,:))),'b-','LineWidth',3),
plot(20*log10(abs(series.H_modCM_2f1_f2(4,:))),'bx-')
plot(20*log10(abs(series.H_modCM_2f1_f2(5,:))),'b.:')
plot(20*log10(abs(series.l_unmodCM_2f1_f2)),'c.','MarkerEdgeColor',[.6 .6 1])
plot(20*log10(abs(series.H_modCM_f2_f1(1,:))),'r.:')
plot(20*log10(abs(series.H_modCM_f2_f1(2,:))),'rx-')
plot(20*log10(abs(series.H_modCM_f2_f1(3,:))),'r-','LineWidth',3),
plot(20*log10(abs(series.H_modCM_f2_f1(4,:))),'rx-')
plot(20*log10(abs(series.H_modCM_f2_f1(5,:))),'r.:')
plot(20*log10(abs(series.l_unmodCM_f2_f1)),'c.','MarkerEdgeColor',[1 .6 .6])
grid on,title modDPCM,set(gca,'FontSize',7); a=a+1;

% Plot modSF level
h_ax(a)=subplot('position',[0.0446    0.41    0.9429    0.18]);
plot(20*log10(abs(series.H_mod_DPf1(1,:))),'g.:','MarkerEdgeColor',[0 .7 0],'Color',[0 .7 0]), hold on
plot(20*log10(abs(series.H_mod_DPf1(2,:))),'gx-','MarkerEdgeColor',[0 .7 0],'Color',[0 .7 0])
plot(20*log10(abs(series.H_mod_DPf1(3,:))),'g-','LineWidth',3,'MarkerEdgeColor',[0 .7 0],'Color',[0 .7 0]),
plot(20*log10(abs(series.H_mod_DPf1(4,:))),'gx-','MarkerEdgeColor',[0 .7 0],'Color',[0 .7 0])
plot(20*log10(abs(series.H_mod_DPf1(5,:))),'g.:','MarkerEdgeColor',[0 .6 0],'Color',[0 .7 0])
plot(20*log10(abs(series.l_unmod_SF_DPf1)),'c.','MarkerEdgeColor',[.6 .8 .6])
plot(20*log10(abs(series.H_mod_DPf2(1,:))),'w.:'), hold on
plot(20*log10(abs(series.H_mod_DPf2(2,:))),'wx-')
plot(20*log10(abs(series.H_mod_DPf2(3,:))),'w-','LineWidth',3),
plot(20*log10(abs(series.H_mod_DPf2(4,:))),'wx-')
plot(20*log10(abs(series.H_mod_DPf2(5,:))),'w.:')
plot(20*log10(abs(series.l_unmod_SF_DPf2)),'c.','MarkerEdgeColor',[.9 .9 .9])
plot(20*log10(abs(series.H_mod_f2(1,:))),'k.:'), hold on
plot(20*log10(abs(series.H_mod_f2(2,:))),'kx-')
plot(20*log10(abs(series.H_mod_f2(3,:))),'k-','LineWidth',3),
plot(20*log10(abs(series.H_mod_f2(4,:))),'kx-')
plot(20*log10(abs(series.H_mod_f2(5,:))),'k.:')
plot(20*log10(abs(series.l_unmod_SF_f2)),'c.','MarkerEdgeColor',[.4 .4 .4])
grid on,title modSF,set(gca,'FontSize',7); a=a+1;

% Plot modCM_f2 level
h_ax(a)=subplot('position',[0.0446   0.222   0.9429    0.18]);
plot(20*log10(abs(series.H_modCM_DPf1(1,:))),'g.:','MarkerEdgeColor',[0 .7 0],'Color',[0 .7 0]), hold on
plot(20*log10(abs(series.H_modCM_DPf1(2,:))),'gx-','MarkerEdgeColor',[0 .7 0],'Color',[0 .7 0])
plot(20*log10(abs(series.H_modCM_DPf1(3,:))),'g-','LineWidth',3,'MarkerEdgeColor',[0 .7 0],'Color',[0 .7 0]),
plot(20*log10(abs(series.H_modCM_DPf1(4,:))),'gx-','MarkerEdgeColor',[0 .7 0],'Color',[0 .7 0])
plot(20*log10(abs(series.H_modCM_DPf1(5,:))),'g.:','MarkerEdgeColor',[0 .6 0],'Color',[0 .7 0])
plot(20*log10(abs(series.l_unmodCM_DPf1)),'c.','MarkerEdgeColor',[.6 .8 .6])
plot(20*log10(abs(series.H_modCM_DPf2(1,:))),'w.:'), hold on
plot(20*log10(abs(series.H_modCM_DPf2(2,:))),'wx-')
plot(20*log10(abs(series.H_modCM_DPf2(3,:))),'w-','LineWidth',3),
plot(20*log10(abs(series.H_modCM_DPf2(4,:))),'wx-')
plot(20*log10(abs(series.H_modCM_DPf2(5,:))),'w.:')
plot(20*log10(abs(series.l_unmodCM_DPf2)),'c.','MarkerEdgeColor',[.9 .9 .9])
plot(20*log10(abs(series.H_modCM_f2(1,:))),'k.:'), hold on
plot(20*log10(abs(series.H_modCM_f2(2,:))),'kx-')
plot(20*log10(abs(series.H_modCM_f2(3,:))),'k-','LineWidth',3),
plot(20*log10(abs(series.H_modCM_f2(4,:))),'kx-')
plot(20*log10(abs(series.H_modCM_f2(5,:))),'k.:')
plot(20*log10(abs(series.l_unmodCM_f2)),'c.','MarkerEdgeColor',[.4 .4 .4])
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
xlim(h_ax(1),[1 len_series])
zoom on
