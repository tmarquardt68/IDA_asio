function plot_time_phase_series_spectral_lines(series,x,current_fig)
% To do
% !!! Phases from 2nd entry onwards becasue often no BT in 1st entry. Make an option in the GUI?
% plot CAP timing (latency re. expected)

len_series = length(series.l_unmod_SF_f2);

a=1;
% Plot modDP level
h_ax(a)=subplot('position',[0.0446    0.811    0.9429    0.18]);
% plot(t,[NaN (angle(series.H_mod_2f2_f1(1,2:end-1)/series.H_mod_2f2_f1(1,2)))/2/pi NaN],'y.:','LineWidth',1), hold on
% plot(t,[NaN (angle(series.H_mod_2f2_f1(2,2:end-1)/series.H_mod_2f2_f1(2,2)))/2/pi NaN],'yx-')
% plot(t,[NaN (angle(series.H_mod_2f2_f1(3,2:end-1)/series.H_mod_2f2_f1(3,2)))/2/pi NaN],'y-','LineWidth',3),
% plot(t,[NaN (angle(series.H_mod_2f2_f1(4,2:end-1)/series.H_mod_2f2_f1(4,2)))/2/pi NaN],'yx-')
% plot(t,[NaN (angle(series.H_mod_2f2_f1(5,2:end-1)/series.H_mod_2f2_f1(5,2)))/2/pi NaN],'y.:','LineWidth',1)
plot(x,(angle(series.l_unmod_2f2_f1/series.l_unmod_2f2_f1(1,2)))/2/pi,'c.','MarkerEdgeColor',[.9 .9 .5])
plot(x,[NaN (angle(series.H_mod_2f1(1,2:end-1)/series.H_mod_2f1(1,2)))/2/pi NaN],'m.:','LineWidth',1), hold on
plot(x,[NaN (angle(series.H_mod_2f1(2,2:end-1)/series.H_mod_2f1(2,2)))/2/pi NaN],'mx-')
plot(x,[NaN (angle(series.H_mod_2f1(3,2:end-1)/series.H_mod_2f1(3,2)))/2/pi NaN],'m-','LineWidth',3),
plot(x,[NaN (angle(series.H_mod_2f1(4,2:end-1)/series.H_mod_2f1(4,2)))/2/pi NaN],'mx-')
plot(x,[NaN (angle(series.H_mod_2f1(5,2:end-1)/series.H_mod_2f1(5,2)))/2/pi NaN],'m.:','LineWidth',1) 
% plot(t,(angle(series.l_unmod_2f1)/series.H_mod_2f2_f1(1,2)))/2/pi,'c.','MarkerEdgeColor',[1 .6 1]) %curently not a filed in $results
plot(x,[NaN (angle(series.H_mod_2f1_f2(1,2:end-1)/series.H_mod_2f1_f2(1,2)))/2/pi NaN],'b.:','LineWidth',1), hold on
plot(x,[NaN (angle(series.H_mod_2f1_f2(2,2:end-1)/series.H_mod_2f1_f2(2,2)))/2/pi NaN],'bx-')
plot(x,[NaN (angle(series.H_mod_2f1_f2(3,2:end-1)/series.H_mod_2f1_f2(3,2)))/2/pi NaN],'b-','LineWidth',3),
plot(x,[NaN (angle(series.H_mod_2f1_f2(4,2:end-1)/series.H_mod_2f1_f2(4,2)))/2/pi NaN],'bx-')
plot(x,[NaN (angle(series.H_mod_2f1_f2(5,2:end-1)/series.H_mod_2f1_f2(5,2)))/2/pi NaN],'b.:','LineWidth',1)
plot(x,(angle(series.l_unmod_2f1_f2/series.l_unmod_2f1_f2(1,2)))/2/pi,'c.','MarkerEdgeColor',[.6 .6 1])
plot(x,[NaN (angle(series.H_mod_f2_f1(1,2:end-1)/series.H_mod_f2_f1(1,2)))/2/pi NaN],'r.:','LineWidth',1)
plot(x,[NaN (angle(series.H_mod_f2_f1(2,2:end-1)/series.H_mod_f2_f1(2,2)))/2/pi NaN],'rx-')
plot(x,[NaN (angle(series.H_mod_f2_f1(3,2:end-1)/series.H_mod_f2_f1(3,2)))/2/pi NaN],'r-','LineWidth',3),
plot(x,[NaN (angle(series.H_mod_f2_f1(4,2:end-1)/series.H_mod_f2_f1(4,2)))/2/pi NaN],'rx-')
plot(x,[NaN (angle(series.H_mod_f2_f1(5,2:end-1)/series.H_mod_f2_f1(5,2)))/2/pi NaN],'r.:','LineWidth',1)
plot(x,(angle(series.l_unmod_f2_f1/series.l_unmod_f2_f1(1,2)))/2/pi,'c.','MarkerEdgeColor',[1 .6 .6])
if isduration(x),xtickformat('mm:ss'),end ,grid on,title modDP,set(gca,'FontSize',7);a=a+1;
% Plot modDP CM
h_ax(a)=subplot('position',[0.0446    0.617    0.9429    0.18]);
% plot(t,[NaN (angle(series.H_modCM_2f2_f1(1,2:end-1)/series.H_modCM_2f2_f1(1,2)))/2/pi NaN],'y.:','LineWidth',1), hold on
% plot(t,[NaN (angle(series.H_modCM_2f2_f1(2,2:end-1)/series.H_modCM_2f2_f1(2,2)))/2/pi NaN],'yx-')
% plot(t,[NaN (angle(series.H_modCM_2f2_f1(3,2:end-1)/series.H_modCM_2f2_f1(3,2)))/2/pi NaN],'y-','LineWidth',3),
% plot(t,[NaN (angle(series.H_modCM_2f2_f1(4,2:end-1)/series.H_modCM_2f2_f1(4,2)))/2/pi NaN],'yx-')
% plot(t,[NaN (angle(series.H_modCM_2f2_f1(5,2:end-1)/series.H_modCM_2f2_f1(5,2)))/2/pi NaN],'y.:','LineWidth',1)
plot(x,(angle(series.l_unmodCM_2f2_f1/series.l_unmodCM_2f2_f1(1,2)))/2/pi,'c.','MarkerEdgeColor',[.9 .9 .5])
plot(x,[NaN (angle(series.H_modCM_2f1(1,2:end-1)/series.H_modCM_2f1(1,2)))/2/pi NaN],'m.:','LineWidth',1), hold on
plot(x,[NaN (angle(series.H_modCM_2f1(2,2:end-1)/series.H_modCM_2f1(2,2)))/2/pi NaN],'mx-')
plot(x,[NaN (angle(series.H_modCM_2f1(3,2:end-1)/series.H_modCM_2f1(3,2)))/2/pi NaN],'m-','LineWidth',3),
plot(x,[NaN (angle(series.H_modCM_2f1(4,2:end-1)/series.H_modCM_2f1(4,2)))/2/pi NaN],'mx-')
plot(x,[NaN (angle(series.H_modCM_2f1(5,2:end-1)/series.H_modCM_2f1(5,2)))/2/pi NaN],'m.:','LineWidth',1)
plot(x,(angle(series.l_unmodCM_2f1/series.l_unmodCM_2f1(1,2)))/2/pi,'c.','MarkerEdgeColor',[1 .6 1])
plot(x,[NaN (angle(series.H_modCM_2f1_f2(1,2:end-1)/series.H_modCM_2f1_f2(1,2)))/2/pi NaN],'b.:','LineWidth',1), hold on
plot(x,[NaN (angle(series.H_modCM_2f1_f2(2,2:end-1)/series.H_modCM_2f1_f2(2,2)))/2/pi NaN],'bx-')
plot(x,[NaN (angle(series.H_modCM_2f1_f2(3,2:end-1)/series.H_modCM_2f1_f2(3,2)))/2/pi NaN],'b-','LineWidth',3),
plot(x,[NaN (angle(series.H_modCM_2f1_f2(4,2:end-1)/series.H_modCM_2f1_f2(4,2)))/2/pi NaN],'bx-')
plot(x,[NaN (angle(series.H_modCM_2f1_f2(5,2:end-1)/series.H_modCM_2f1_f2(5,2)))/2/pi NaN],'b.:','LineWidth',1)
plot(x,(angle(series.l_unmodCM_2f1_f2/series.l_unmodCM_2f1_f2(1,2)))/2/pi,'c.','MarkerEdgeColor',[.6 .6 1])
plot(x,[NaN (angle(series.H_modCM_f2_f1(1,2:end-1)/series.H_modCM_f2_f1(1,2)))/2/pi NaN],'r.:','LineWidth',1)
plot(x,[NaN (angle(series.H_modCM_f2_f1(2,2:end-1)/series.H_modCM_f2_f1(2,2)))/2/pi NaN],'rx-')
plot(x,[NaN (angle(series.H_modCM_f2_f1(3,2:end-1)/series.H_modCM_f2_f1(3,2)))/2/pi NaN],'r-','LineWidth',3),
plot(x,[NaN (angle(series.H_modCM_f2_f1(4,2:end-1)/series.H_modCM_f2_f1(4,2)))/2/pi NaN],'rx-')
plot(x,[NaN (angle(series.H_modCM_f2_f1(5,2:end-1)/series.H_modCM_f2_f1(5,2)))/2/pi NaN],'r.:','LineWidth',1)
plot(x,(angle(series.l_unmodCM_f2_f1/series.l_unmodCM_f2_f1(1,2)))/2/pi,'c.','MarkerEdgeColor',[1 .6 .6])
if isduration(x),xtickformat('mm:ss'),end ,grid on,title modDPCM,set(gca,'FontSize',7); a=a+1;

% Plot modSF level
h_ax(a)=subplot('position',[0.0446    0.423    0.9429    0.18]);
if isfield(series,'H_mod_DPf1')
    plot(x,[NaN (angle(series.H_mod_DPf1(1,2:end-1)/series.H_mod_DPf1(1,2)))/2/pi NaN],'g.:','MarkerEdgeColor',[0 .7 0],'Color',[0 .7 0]), hold on
    plot(x,[NaN (angle(series.H_mod_DPf1(2,2:end-1)/series.H_mod_DPf1(2,2)))/2/pi NaN],'gx-','MarkerEdgeColor',[0 .7 0],'Color',[0 .7 0])
    plot(x,[NaN (angle(series.H_mod_DPf1(3,2:end-1)/series.H_mod_DPf1(3,2)))/2/pi NaN],'g-','LineWidth',3,'MarkerEdgeColor',[0 .7 0],'Color',[0 .7 0]),
    plot(x,[NaN (angle(series.H_mod_DPf1(4,2:end-1)/series.H_mod_DPf1(4,2)))/2/pi NaN],'gx-','MarkerEdgeColor',[0 .7 0],'Color',[0 .7 0])
    plot(x,[NaN (angle(series.H_mod_DPf1(5,2:end-1)/series.H_mod_DPf1(5,2)))/2/pi NaN],'g.:','MarkerEdgeColor',[0 .6 0],'Color',[0 .7 0])
    plot(x,(angle(series.l_unmod_SF_DPf1/series.l_unmod_SF_DPf1(1,2)))/2/pi,'c.','MarkerEdgeColor',[.6 .8 .6])
    plot(x,[NaN (angle(series.H_mod_DPf2(1,2:end-1)/series.H_mod_DPf2(1,2)))/2/pi NaN],'w.:','LineWidth',1), hold on
    plot(x,[NaN (angle(series.H_mod_DPf2(2,2:end-1)/series.H_mod_DPf2(2,2)))/2/pi NaN],'wx-')
    plot(x,[NaN (angle(series.H_mod_DPf2(3,2:end-1)/series.H_mod_DPf2(3,2)))/2/pi NaN],'w-','LineWidth',3),
    plot(x,[NaN (angle(series.H_mod_DPf2(4,2:end-1)/series.H_mod_DPf2(4,2)))/2/pi NaN],'wx-')
    plot(x,[NaN (angle(series.H_mod_DPf2(5,2:end-1)/series.H_mod_DPf2(5,2)))/2/pi NaN],'w.:','LineWidth',1)
    plot(x,(angle(series.l_unmod_SF_DPf2/series.l_unmod_SF_DPf2(1,2)))/2/pi,'c.','MarkerEdgeColor',[.9 .9 .9])
    plot(x,[NaN (angle(series.H_mod_f2(1,2:end-1)/series.H_mod_f2(1,2)))/2/pi NaN],'k.:','LineWidth',1), hold on
    plot(x,[NaN (angle(series.H_mod_f2(2,2:end-1)/series.H_mod_f2(2,2)))/2/pi NaN],'kx-')
    plot(x,[NaN (angle(series.H_mod_f2(3,2:end-1)/series.H_mod_f2(3,2)))/2/pi NaN],'k-','LineWidth',3),
    plot(x,[NaN (angle(series.H_mod_f2(4,2:end-1)/series.H_mod_f2(4,2)))/2/pi NaN],'kx-')
    plot(x,[NaN (angle(series.H_mod_f2(5,2:end-1)/series.H_mod_f2(5,2)))/2/pi NaN],'k.:','LineWidth',1)
    plot(x,(angle(series.l_unmod_SF_f2/series.l_unmod_SF_f2(1,2)))/2/pi,'c.','MarkerEdgeColor',[.4 .4 .4])
    if isduration(x),xtickformat('mm:ss'),end ,grid on,title modSF,set(gca,'FontSize',7);a=a+1;
end

% Plot modCM_f2 level
h_ax(a)=subplot('position',[0.0446   0.228   0.9429    0.18]);
plot(x,[NaN (angle(series.H_modCM_DPf1(1,2:end-1)/series.H_modCM_DPf1(1,2)))/2/pi NaN],'g.:','MarkerEdgeColor',[0 .7 0],'Color',[0 .7 0]), hold on
plot(x,[NaN (angle(series.H_modCM_DPf1(2,2:end-1)/series.H_modCM_DPf1(2,2)))/2/pi NaN],'gx-','MarkerEdgeColor',[0 .7 0],'Color',[0 .7 0])
plot(x,[NaN (angle(series.H_modCM_DPf1(3,2:end-1)/series.H_modCM_DPf1(3,2)))/2/pi NaN],'g-','LineWidth',3,'MarkerEdgeColor',[0 .7 0],'Color',[0 .7 0]),
plot(x,[NaN (angle(series.H_modCM_DPf1(4,2:end-1)/series.H_modCM_DPf1(4,2)))/2/pi NaN],'gx-','MarkerEdgeColor',[0 .7 0],'Color',[0 .7 0])
plot(x,[NaN (angle(series.H_modCM_DPf1(5,2:end-1)/series.H_modCM_DPf1(5,2)))/2/pi NaN],'g.:','MarkerEdgeColor',[0 .6 0],'Color',[0 .7 0])
plot(x,(angle(series.l_unmodCM_DPf1/series.l_unmodCM_DPf1(1,2)))/2/pi,'c.','MarkerEdgeColor',[.6 .8 .6])
plot(x,[NaN (angle(series.H_modCM_DPf2(1,2:end-1)/series.H_modCM_DPf2(1,2)))/2/pi NaN],'w.:','LineWidth',1), hold on
plot(x,[NaN (angle(series.H_modCM_DPf2(2,2:end-1)/series.H_modCM_DPf2(2,2)))/2/pi NaN],'wx-')
plot(x,[NaN (angle(series.H_modCM_DPf2(3,2:end-1)/series.H_modCM_DPf2(3,2)))/2/pi NaN],'w-','LineWidth',3),
plot(x,[NaN (angle(series.H_modCM_DPf2(4,2:end-1)/series.H_modCM_DPf2(4,2)))/2/pi NaN],'wx-')
plot(x,[NaN (angle(series.H_modCM_DPf2(5,2:end-1)/series.H_modCM_DPf2(5,2)))/2/pi NaN],'w.:','LineWidth',1)
plot(x,(angle(series.l_unmodCM_DPf2/series.l_unmodCM_DPf2(1,2)))/2/pi,'c.','MarkerEdgeColor',[.9 .9 .9])
plot(x,[NaN (angle(series.H_modCM_f2(1,2:end-1)/series.H_modCM_f2(1,2)))/2/pi NaN],'k.:','LineWidth',1), hold on
plot(x,[NaN (angle(series.H_modCM_f2(2,2:end-1)/series.H_modCM_f2(2,2)))/2/pi NaN],'kx-')
plot(x,[NaN (angle(series.H_modCM_f2(3,2:end-1)/series.H_modCM_f2(3,2)))/2/pi NaN],'k-','LineWidth',3),
plot(x,[NaN (angle(series.H_modCM_f2(4,2:end-1)/series.H_modCM_f2(4,2)))/2/pi NaN],'kx-')
plot(x,[NaN (angle(series.H_modCM_f2(5,2:end-1)/series.H_modCM_f2(5,2)))/2/pi NaN],'k.:','LineWidth',1)
plot(x,(angle(series.l_unmodCM_f2/series.l_unmodCM_f2(1,2)))/2/pi,'c.','MarkerEdgeColor',[.4 .4 .4])
if isduration(x),xtickformat('mm:ss'),end ,grid on,title modCMf2,set(gca,'FontSize',7);a=a+1;

% Plot CAPs
title CAPs,set(gca,'FontSize',7,'Tag','axes_CAP')
if isfield(series,'CAP_ampl')
    h_ax(a)=subplot('position',[0.0446   0.033   0.9429    0.18]); 
    CAP_amp_series = reshape([series.CAP_ampl],8,len_series)';
    plot(x,CAP_amp_series(1:end,:)) 
    if isduration(x),xtickformat('mm:ss'),end ,grid on,title CAPs,set(gca,'FontSize',7,'Tag','axes_CAP')
end

xlabel('time re. hypoxia onset [MM:SS]'),zoom on
linkaxes(h_ax, 'x')
set(h_ax,'Color',[.8 .8 .8],'FontSize',7)
xlim(h_ax(1),[x(1) x(end)])
zoom on

setappdata(current_fig,'h_axes',h_ax)