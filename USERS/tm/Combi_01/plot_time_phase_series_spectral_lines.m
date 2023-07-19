function plot_time_phase_series_spectral_lines(series)

a=1;
% Plot modDP level
h_ax(a)=subplot('position',[0.0446    0.79    0.9429    0.18]);
plot((angle(series.H_mod_2f2_f1(1,2:end)/series.H_mod_2f2_f1(1,2)))/2/pi,'y.:'), hold on
plot((angle(series.H_mod_2f2_f1(2,2:end)/series.H_mod_2f2_f1(2,2)))/2/pi,'yx-')
plot((angle(series.H_mod_2f2_f1(3,2:end)/series.H_mod_2f2_f1(3,2)))/2/pi,'y-','LineWidth',3),
plot((angle(series.H_mod_2f2_f1(4,2:end)/series.H_mod_2f2_f1(4,2)))/2/pi,'yx-')
plot((angle(series.H_mod_2f2_f1(5,2:end)/series.H_mod_2f2_f1(5,2)))/2/pi,'y.:')
plot((angle(series.l_unmod_2f2_f1/series.l_unmod_2f2_f1(1,2)))/2/pi,'c.','MarkerEdgeColor',[.9 .9 .5])
plot((angle(series.H_mod_2f1(1,2:end)/series.H_mod_2f1(1,2)))/2/pi,'m.:'), hold on
plot((angle(series.H_mod_2f1(2,2:end)/series.H_mod_2f1(2,2)))/2/pi,'mx-')
plot((angle(series.H_mod_2f1(3,2:end)/series.H_mod_2f1(3,2)))/2/pi,'m-','LineWidth',3),
plot((angle(series.H_mod_2f1(4,2:end)/series.H_mod_2f1(4,2)))/2/pi,'mx-')
plot((angle(series.H_mod_2f1(5,2:end)/series.H_mod_2f1(5,2)))/2/pi,'m.:') 
% plot((angle(series.l_unmod_2f1)/series.H_mod_2f2_f1(1,2)))/2/pi,'c.','MarkerEdgeColor',[1 .6 1]) %curently not a filed in $results
plot((angle(series.H_mod_2f1_f2(1,2:end)/series.H_mod_2f1_f2(1,2)))/2/pi,'b.:'), hold on
plot((angle(series.H_mod_2f1_f2(2,2:end)/series.H_mod_2f1_f2(2,2)))/2/pi,'bx-')
plot((angle(series.H_mod_2f1_f2(3,2:end)/series.H_mod_2f1_f2(3,2)))/2/pi,'b-','LineWidth',3),
plot((angle(series.H_mod_2f1_f2(4,2:end)/series.H_mod_2f1_f2(4,2)))/2/pi,'bx-')
plot((angle(series.H_mod_2f1_f2(5,2:end)/series.H_mod_2f1_f2(5,2)))/2/pi,'b.:')
plot((angle(series.l_unmod_2f1_f2/series.l_unmod_2f1_f2(1,2)))/2/pi,'c.','MarkerEdgeColor',[.6 .6 1])
plot((angle(series.H_mod_f2_f1(1,2:end)/series.H_mod_f2_f1(1,2)))/2/pi,'r.:')
plot((angle(series.H_mod_f2_f1(2,2:end)/series.H_mod_f2_f1(2,2)))/2/pi,'rx-')
plot((angle(series.H_mod_f2_f1(3,2:end)/series.H_mod_f2_f1(3,2)))/2/pi,'r-','LineWidth',3),
plot((angle(series.H_mod_f2_f1(4,2:end)/series.H_mod_f2_f1(4,2)))/2/pi,'rx-')
plot((angle(series.H_mod_f2_f1(5,2:end)/series.H_mod_f2_f1(5,2)))/2/pi,'r.:')
plot((angle(series.l_unmod_f2_f1/series.l_unmod_f2_f1(1,2)))/2/pi,'c.','MarkerEdgeColor',[1 .6 .6])
grid on,title modDP,set(gca,'FontSize',7); a=a+1;

% Plot modDP CM
h_ax(a)=subplot('position',[0.0446    0.60    0.9429    0.18]);
plot((angle(series.H_modCM_2f2_f1(1,2:end)/series.H_modCM_2f2_f1(1,2)))/2/pi,'y.:'), hold on
plot((angle(series.H_modCM_2f2_f1(2,2:end)/series.H_modCM_2f2_f1(2,2)))/2/pi,'yx-')
plot((angle(series.H_modCM_2f2_f1(3,2:end)/series.H_modCM_2f2_f1(3,2)))/2/pi,'y-','LineWidth',3),
plot((angle(series.H_modCM_2f2_f1(4,2:end)/series.H_modCM_2f2_f1(4,2)))/2/pi,'yx-')
plot((angle(series.H_modCM_2f2_f1(5,2:end)/series.H_modCM_2f2_f1(5,2)))/2/pi,'y.:')
plot((angle(series.l_unmodCM_2f2_f1/series.l_unmodCM_2f2_f1(1,2)))/2/pi,'c.','MarkerEdgeColor',[.9 .9 .5])
plot((angle(series.H_modCM_2f1(1,2:end)/series.H_modCM_2f1(1,2)))/2/pi,'m.:'), hold on
plot((angle(series.H_modCM_2f1(2,2:end)/series.H_modCM_2f1(2,2)))/2/pi,'mx-')
plot((angle(series.H_modCM_2f1(3,2:end)/series.H_modCM_2f1(3,2)))/2/pi,'m-','LineWidth',3),
plot((angle(series.H_modCM_2f1(4,2:end)/series.H_modCM_2f1(4,2)))/2/pi,'mx-')
plot((angle(series.H_modCM_2f1(5,2:end)/series.H_modCM_2f1(5,2)))/2/pi,'m.:')
plot((angle(series.l_unmodCM_2f1/series.l_unmodCM_2f1(1,2)))/2/pi,'c.','MarkerEdgeColor',[1 .6 1])
plot((angle(series.H_modCM_2f1_f2(1,2:end)/series.H_modCM_2f1_f2(1,2)))/2/pi,'b.:'), hold on
plot((angle(series.H_modCM_2f1_f2(2,2:end)/series.H_modCM_2f1_f2(2,2)))/2/pi,'bx-')
plot((angle(series.H_modCM_2f1_f2(3,2:end)/series.H_modCM_2f1_f2(3,2)))/2/pi,'b-','LineWidth',3),
plot((angle(series.H_modCM_2f1_f2(4,2:end)/series.H_modCM_2f1_f2(4,2)))/2/pi,'bx-')
plot((angle(series.H_modCM_2f1_f2(5,2:end)/series.H_modCM_2f1_f2(5,2)))/2/pi,'b.:')
plot((angle(series.l_unmodCM_2f1_f2/series.l_unmodCM_2f1_f2(1,2)))/2/pi,'c.','MarkerEdgeColor',[.6 .6 1])
plot((angle(series.H_modCM_f2_f1(1,2:end)/series.H_modCM_f2_f1(1,2)))/2/pi,'r.:')
plot((angle(series.H_modCM_f2_f1(2,2:end)/series.H_modCM_f2_f1(2,2)))/2/pi,'rx-')
plot((angle(series.H_modCM_f2_f1(3,2:end)/series.H_modCM_f2_f1(3,2)))/2/pi,'r-','LineWidth',3),
plot((angle(series.H_modCM_f2_f1(4,2:end)/series.H_modCM_f2_f1(4,2)))/2/pi,'rx-')
plot((angle(series.H_modCM_f2_f1(5,2:end)/series.H_modCM_f2_f1(5,2)))/2/pi,'r.:')
plot((angle(series.l_unmodCM_f2_f1/series.l_unmodCM_f2_f1(1,2)))/2/pi,'c.','MarkerEdgeColor',[1 .6 .6])
grid on,title modDPCM,set(gca,'FontSize',7); a=a+1;

% Plot modSF level
h_ax(a)=subplot('position',[0.0446    0.41    0.9429    0.18]);
plot((angle(series.H_mod_DPf1(1,2:end)/series.H_mod_DPf1(1,2)))/2/pi,'g.:','MarkerEdgeColor',[0 .7 0],'Color',[0 .7 0]), hold on
plot((angle(series.H_mod_DPf1(2,2:end)/series.H_mod_DPf1(2,2)))/2/pi,'gx-','MarkerEdgeColor',[0 .7 0],'Color',[0 .7 0])
plot((angle(series.H_mod_DPf1(3,2:end)/series.H_mod_DPf1(3,2)))/2/pi,'g-','LineWidth',3,'MarkerEdgeColor',[0 .7 0],'Color',[0 .7 0]),
plot((angle(series.H_mod_DPf1(4,2:end)/series.H_mod_DPf1(4,2)))/2/pi,'gx-','MarkerEdgeColor',[0 .7 0],'Color',[0 .7 0])
plot((angle(series.H_mod_DPf1(5,2:end)/series.H_mod_DPf1(5,2)))/2/pi,'g.:','MarkerEdgeColor',[0 .6 0],'Color',[0 .7 0])
plot((angle(series.l_unmod_SF_DPf1/series.l_unmod_SF_DPf1(1,2)))/2/pi,'c.','MarkerEdgeColor',[.6 .8 .6])
plot((angle(series.H_mod_DPf2(1,2:end)/series.H_mod_DPf2(1,2)))/2/pi,'w.:'), hold on
plot((angle(series.H_mod_DPf2(2,2:end)/series.H_mod_DPf2(2,2)))/2/pi,'wx-')
plot((angle(series.H_mod_DPf2(3,2:end)/series.H_mod_DPf2(3,2)))/2/pi,'w-','LineWidth',3),
plot((angle(series.H_mod_DPf2(4,2:end)/series.H_mod_DPf2(4,2)))/2/pi,'wx-')
plot((angle(series.H_mod_DPf2(5,2:end)/series.H_mod_DPf2(5,2)))/2/pi,'w.:')
plot((angle(series.l_unmod_SF_DPf2/series.l_unmod_SF_DPf2(1,2)))/2/pi,'c.','MarkerEdgeColor',[.9 .9 .9])
plot((angle(series.H_mod_f2(1,2:end)/series.H_mod_f2(1,2)))/2/pi,'k.:'), hold on
plot((angle(series.H_mod_f2(2,2:end)/series.H_mod_f2(2,2)))/2/pi,'kx-')
plot((angle(series.H_mod_f2(3,2:end)/series.H_mod_f2(3,2)))/2/pi,'k-','LineWidth',3),
plot((angle(series.H_mod_f2(4,2:end)/series.H_mod_f2(4,2)))/2/pi,'kx-')
plot((angle(series.H_mod_f2(5,2:end)/series.H_mod_f2(5,2)))/2/pi,'k.:')
plot((angle(series.l_unmod_SF_f2/series.l_unmod_SF_f2(1,2)))/2/pi,'c.','MarkerEdgeColor',[.4 .4 .4])
grid on,title modSF,set(gca,'FontSize',7); a=a+1;

% Plot modCM_f2 level
h_ax(a)=subplot('position',[0.0446   0.222   0.9429    0.18]);
plot((angle(series.H_modCM_DPf1(1,2:end)/series.H_modCM_DPf1(1,2)))/2/pi,'g.:','MarkerEdgeColor',[0 .7 0],'Color',[0 .7 0]), hold on
plot((angle(series.H_modCM_DPf1(2,2:end)/series.H_modCM_DPf1(2,2)))/2/pi,'gx-','MarkerEdgeColor',[0 .7 0],'Color',[0 .7 0])
plot((angle(series.H_modCM_DPf1(3,2:end)/series.H_modCM_DPf1(3,2)))/2/pi,'g-','LineWidth',3,'MarkerEdgeColor',[0 .7 0],'Color',[0 .7 0]),
plot((angle(series.H_modCM_DPf1(4,2:end)/series.H_modCM_DPf1(4,2)))/2/pi,'gx-','MarkerEdgeColor',[0 .7 0],'Color',[0 .7 0])
plot((angle(series.H_modCM_DPf1(5,2:end)/series.H_modCM_DPf1(5,2)))/2/pi,'g.:','MarkerEdgeColor',[0 .6 0],'Color',[0 .7 0])
plot((angle(series.l_unmodCM_DPf1/series.l_unmodCM_DPf1(1,2)))/2/pi,'c.','MarkerEdgeColor',[.6 .8 .6])
plot((angle(series.H_modCM_DPf2(1,2:end)/series.H_modCM_DPf2(1,2)))/2/pi,'w.:'), hold on
plot((angle(series.H_modCM_DPf2(2,2:end)/series.H_modCM_DPf2(2,2)))/2/pi,'wx-')
plot((angle(series.H_modCM_DPf2(3,2:end)/series.H_modCM_DPf2(3,2)))/2/pi,'w-','LineWidth',3),
plot((angle(series.H_modCM_DPf2(4,2:end)/series.H_modCM_DPf2(4,2)))/2/pi,'wx-')
plot((angle(series.H_modCM_DPf2(5,2:end)/series.H_modCM_DPf2(5,2)))/2/pi,'w.:')
plot((angle(series.l_unmodCM_DPf2/series.l_unmodCM_DPf2(1,2)))/2/pi,'c.','MarkerEdgeColor',[.9 .9 .9])
plot((angle(series.H_modCM_f2(1,2:end)/series.H_modCM_f2(1,2)))/2/pi,'k.:'), hold on
plot((angle(series.H_modCM_f2(2,2:end)/series.H_modCM_f2(2,2)))/2/pi,'kx-')
plot((angle(series.H_modCM_f2(3,2:end)/series.H_modCM_f2(3,2)))/2/pi,'k-','LineWidth',3),
plot((angle(series.H_modCM_f2(4,2:end)/series.H_modCM_f2(4,2)))/2/pi,'kx-')
plot((angle(series.H_modCM_f2(5,2:end)/series.H_modCM_f2(5,2)))/2/pi,'k.:')
plot((angle(series.l_unmodCM_f2/series.l_unmodCM_f2(1,2)))/2/pi,'c.','MarkerEdgeColor',[.4 .4 .4])
grid on,title modCMf2,set(gca,'FontSize',7); a=a+1;

% Plot CAPs
len_series = length(series.l_unmod_SF_f2);
h_ax(a)=subplot('position',[0.0446   0.03   0.9429    0.18]); 
%plot(reshape([series.CAP_ampl],8,len_series)'+repmat([0 1 2 3 4 5 6 7]* 0.05,len_series,1))
CAP_amp_series = reshape([series.CAP_ampl],8,len_series)';
plot(CAP_amp_series(2:end,:))
grid on,title CAPs,set(gca,'FontSize',7)

xlabel('no of 6s-presentions'),zoom on
linkaxes(h_ax, 'x')
set(h_ax,'Color',[.8 .8 .8])
xlim(h_ax(1),[1 len_series])
zoom on
