function plot_time_op_series(series,x,current_fig)

len_series = length(series.l_unmod_SF_f2);

a=1;
% Plot modDP level
h_ax(a)=subplot('position',[0.0446    0.811    0.9429    0.18]);

% plot(x,20*log10(abs(series.H_mod_2f1(1,:))+abs(series.H_mod_2f1(5,:))),'m.:','LineWidth',1)
% plot(x,20*log10(abs(series.H_mod_2f1(2,:))+abs(series.H_mod_2f1(4,:))),'mx-')
plot(x,20*log10(abs(series.H_mod_2f1(3,:))),'m-','LineWidth',1)
hold on
plot(x,20*log10(abs(series.l_unmod_2f1)),'m.','MarkerEdgeColor',[1 .6 1]) 
plot(x,20*log10((abs(series.H_mod_2f1(2,:))+abs(series.H_mod_2f1(4,:)))./...
    (abs(series.H_mod_2f1(1,:))+abs(series.H_mod_2f1(5,:)))),'m-','LineWidth',3) 

% plot(x,20*log10(abs(series.H_mod_2f1_f2(1,:))+abs(series.H_mod_2f1_f2(5,:))),'b.:','LineWidth',1), hold on
% plot(x,20*log10(abs(series.H_mod_2f1_f2(2,:))+abs(series.H_mod_2f1_f2(4,:))),'bx-')
plot(x,20*log10(abs(series.H_mod_2f1_f2(3,:))),'b-','LineWidth',1)
hold on
plot(x,20*log10(abs(series.l_unmod_2f1_f2)),'b.','MarkerEdgeColor',[.6 .6 1])
plot(x,20*log10((abs(series.H_mod_2f1_f2(2,:))+abs(series.H_mod_2f1_f2(4,:)))./...
    (abs(series.H_mod_2f1_f2(1,:))+abs(series.H_mod_2f1_f2(5,:)))),'b-','LineWidth',3) 

% plot(x,20*log10(abs(series.H_mod_f2_f1(1,:))+abs(series.H_mod_f2_f1(5,:))),'r.:','LineWidth',1), hold on
% plot(x,20*log10(abs(series.H_mod_f2_f1(2,:))+abs(series.H_mod_f2_f1(4,:))),'rx-')
plot(x,20*log10(abs(series.H_mod_f2_f1(3,:))),'r-','LineWidth',1)
plot(x,20*log10(abs(series.l_unmod_f2_f1)),'r.','MarkerEdgeColor',[1 .6 .6])
plot(x,20*log10((abs(series.H_mod_f2_f1(2,:))+abs(series.H_mod_f2_f1(4,:)))./...
    (abs(series.H_mod_f2_f1(1,:))+abs(series.H_mod_f2_f1(5,:)))),'r-','LineWidth',3) 

plot(x,20*log10(abs(series.l_unmod_2f1_f2)./abs(series.l_unmod_f2_f1)),'k-','LineWidth',3)

plot([x(1) x(end)],[0 0],'k--','LineWidth',1),
if isduration(x),xtickformat('mm:ss'),end ,grid on,title modDP,set(gca,'FontSize',7);a=a+1;y=ylim;ylim([-30,y(2)])

% Plot modDP CM
h_ax(a)=subplot('position',[0.0446    0.617    0.9429    0.18]);

% plot(x,20*log10(abs(series.H_modCM_2f1(1,:))+abs(series.H_modCM_2f1(5,:))),'m.:','LineWidth',1), hold on
% plot(x,20*log10(abs(series.H_modCM_2f1(2,:))+abs(series.H_modCM_2f1(4,:))),'mx-')
plot(x,20*log10(abs(series.H_modCM_2f1(3,:))),'m-','LineWidth',1)
hold on
plot(x,20*log10(abs(series.l_unmodCM_2f1)),'m.','MarkerEdgeColor',[1 .6 1]) 
plot(x,20*log10((abs(series.H_modCM_2f1(2,:))+abs(series.H_modCM_2f1(4,:)))./...
    (abs(series.H_modCM_2f1(1,:))+abs(series.H_modCM_2f1(5,:)))),'m-','LineWidth',3) 

% plot(x,20*log10(abs(series.H_modCM_2f1_f2(1,:))+abs(series.H_modCM_2f1_f2(5,:))),'b.:','LineWidth',1), hold on
% plot(x,20*log10(abs(series.H_modCM_2f1_f2(2,:))+abs(series.H_modCM_2f1_f2(4,:))),'bx-')
plot(x,20*log10(abs(series.H_modCM_2f1_f2(3,:))),'b-','LineWidth',1)
hold on
plot(x,20*log10(abs(series.l_unmodCM_2f1_f2)),'b.','MarkerEdgeColor',[.6 .6 1])
plot(x,20*log10((abs(series.H_modCM_2f1_f2(2,:))+abs(series.H_modCM_2f1_f2(4,:)))./...
    (abs(series.H_modCM_2f1_f2(1,:))+abs(series.H_modCM_2f1_f2(5,:)))),'b-','LineWidth',3) 

% plot(x,20*log10(abs(series.H_modCM_f2_f1(1,:))+abs(series.H_modCM_f2_f1(5,:))),'r.:','LineWidth',1), hold on
% plot(x,20*log10(abs(series.H_modCM_f2_f1(2,:))+abs(series.H_modCM_f2_f1(4,:))),'rx-')
plot(x,20*log10(abs(series.H_modCM_f2_f1(3,:))),'r-','LineWidth',1)
hold on
plot(x,20*log10(abs(series.l_unmodCM_f2_f1)),'r.','MarkerEdgeColor',[1 .6 .6])
plot(x,20*log10((abs(series.H_modCM_f2_f1(2,:))+abs(series.H_modCM_f2_f1(4,:)))./...
    (abs(series.H_modCM_f2_f1(1,:))+abs(series.H_modCM_f2_f1(5,:)))),'r-','LineWidth',3)

plot(x,20*log10(abs(series.l_unmodCM_2f1_f2)./abs(series.l_unmodCM_f2_f1)),'k-','LineWidth',3)


plot([x(1) x(end)],[0 0],'k--','LineWidth',1),
if isduration(x),xtickformat('mm:ss'),end ,grid on,title modDPCM,set(gca,'FontSize',7); a=a+1;

% Plot modSF levels
h_ax(a)=subplot('position',[0.0446    0.423    0.9429    0.18]);
if isfield(series,'H_mod_DPf1')
    plot(x,20*log10(abs(series.H_mod_DPf1(3,:))),'g-','LineWidth',1)
    hold on
    plot(x,20*log10(abs(series.l_unmod_SF_DPf1)),'g.','MarkerEdgeColor',[.6 .8 .6])
    plot(x,20*log10((abs(series.H_mod_DPf1(2,:))+abs(series.H_mod_DPf1(4,:)))./...
        (abs(series.H_mod_DPf1(1,:))+abs(series.H_mod_DPf1(5,:)))),'g-','LineWidth',3)

    plot(x,20*log10(abs(series.H_mod_DPf2(3,:))),'w-','LineWidth',1)
    hold on
    plot(x,20*log10(abs(series.l_unmod_SF_DPf2)),'w.','MarkerEdgeColor',[.9 .9 .9])
    plot(x,20*log10((abs(series.H_mod_DPf2(2,:))+abs(series.H_mod_DPf2(4,:)))./...
        (abs(series.H_mod_DPf2(1,:))+abs(series.H_mod_DPf2(5,:)))),'w-','LineWidth',3)

    plot(x,20*log10(abs(series.H_mod_f2(3,:))),'k-','LineWidth',1)
    hold on
    plot(x,20*log10((abs(series.H_mod_f2(2,:))+abs(series.H_mod_f2(4,:)))./...
        (abs(series.H_mod_f2(1,:))+abs(series.H_mod_f2(5,:)))),'k-','LineWidth',3)
    plot(x,20*log10(abs(series.l_unmod_SF_f2)),'k.','MarkerEdgeColor',[.4 .4 .4])

    plot([x(1) x(end)],[0 0],'k--','LineWidth',1),
    if isduration(x),xtickformat('mm:ss'),end ,grid on,title modSF,set(gca,'FontSize',7);a=a+1;
end


% Plot modCM_f2 level
h_ax(a)=subplot('position',[0.0446   0.228   0.9429    0.18]);

plot(x,20*log10(abs(series.H_modCM_DPf1(3,:))),'g-','LineWidth',1)
hold on
plot(x,20*log10((abs(series.H_modCM_DPf1(2,:))+abs(series.H_modCM_DPf1(4,:)))./...
    (abs(series.H_modCM_DPf1(1,:))+abs(series.H_modCM_DPf1(5,:)))),'g-','LineWidth',3)
plot(x,20*log10(abs(series.l_unmodCM_DPf1)),'g.','MarkerEdgeColor',[.6 .8 .6])

plot(x,20*log10(abs(series.H_modCM_DPf2(3,:))),'w-','LineWidth',1)
hold on
plot(x,20*log10((abs(series.H_modCM_DPf2(2,:))+abs(series.H_modCM_DPf2(4,:)))./...
    (abs(series.H_modCM_DPf2(1,:))+abs(series.H_modCM_DPf2(5,:)))),'w-','LineWidth',3)
plot(x,20*log10(abs(series.l_unmodCM_DPf2)),'w.','MarkerEdgeColor',[.9 .9 .9])

plot(x,20*log10(abs(series.H_modCM_f2(3,:))),'k-','LineWidth',1)
hold on
plot(x,20*log10((abs(series.H_modCM_f2(2,:))+abs(series.H_modCM_f2(4,:)))./...
    (abs(series.H_modCM_f2(1,:))+abs(series.H_modCM_f2(5,:)))),'k-','LineWidth',3)
plot(x,20*log10(abs(series.l_unmodCM_f2)),'k.','MarkerEdgeColor',[.4 .4 .4])

plot([x(1) x(end)],[0 0],'k--','LineWidth',1),
if isduration(x),xtickformat('mm:ss'),end ,grid on,title modCMf2,set(gca,'FontSize',7); a=a+1;

% Plot CAPs
%plot(t,reshape([series.CAP_ampl],8,len_series)'+repmat([0 1 2 3 4 5 6 7]* 0.05,len_series,1))
if isfield(series,'CAP_ampl')
    h_ax(a)=subplot('position',[0.0446   0.033   0.9429    0.18]); 
    CAP_amp_series = reshape([series.CAP_ampl],8,len_series)';
    plot(x,CAP_amp_series(1:end,:)) 
    if isduration(x),xtickformat('mm:ss'),end ,grid on,title CAPs,set(gca,'FontSize',7,'Tag','axes_CAP')
end

xlabel('time re. hypoxia onset [MM:SS]'),zoom on
linkaxes(h_ax, 'x')
set(h_ax,'Color',[.8 .8 .8])
xlim(h_ax(1),[x(1) x(end)])
zoom on

setappdata(current_fig,'h_axes',h_ax)