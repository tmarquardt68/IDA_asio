function plot_time_series_modulation(series)

course_len = series(1).mod_2f1_f2_course;

a=1;
% Plot 2f2-f1 modulation
h_ax(a)=subplot('position',[0.0446    0.80    0.9429    0.18]);
surf([series(:).mod_2f1_f2_course],'EdgeAlpha',0)
set(h_ax(a),'View',[0 -90]),pos=get(h_ax(a),'Position');set(h_ax(a),'Position',[pos(1) ,pos(2),.88,pos(4)]), 
h = colorbar;pos=get(h,'Position');set(h,'Position',[.94 ,pos(2),pos(3)/4,pos(4)])
grid on,title 2f1-f1,set(gca,'FontSize',7); a=a+1;

% Plot modDP CM
h_ax(a)=subplot('position',[0.0446    0.61    0.9429    0.18]);
surf([series(:).modCM_2f1_f2_course],'EdgeAlpha',0)
set(h_ax(a),'View',[0 -90]),pos=get(h_ax(a),'Position');set(h_ax(a),'Position',[pos(1) ,pos(2),.88,pos(4)]), 
h = colorbar;pos=get(h,'Position');set(h,'Position',[.94 ,pos(2),pos(3)/4,pos(4)])
grid on,title CM_2f1-f2,set(gca,'FontSize',7); a=a+1;

% Plot modSF level
h_ax(a)=subplot('position',[0.0446    0.42    0.9429    0.18]);
surf([series(:).mod_f2_f1_course],'EdgeAlpha',0)
set(h_ax(a),'View',[0 -90]),pos=get(h_ax(a),'Position');set(h_ax(a),'Position',[pos(1) ,pos(2),.88,pos(4)]), 
h = colorbar;pos=get(h,'Position');set(h,'Position',[.94 ,pos(2),pos(3)/4,pos(4)])
grid on,title f2-f1,set(gca,'FontSize',7); a=a+1;

% Plot modCM_f2 level
h_ax(a)=subplot('position',[0.0446   0.23   0.9429    0.18]);
surf([series(:).modCM_f2_f1_course],'EdgeAlpha',0)
set(h_ax(a),'View',[0 -90]),pos=get(h_ax(a),'Position');set(h_ax(a),'Position',[pos(1) ,pos(2),.88,pos(4)]), 
h = colorbar;pos=get(h,'Position');set(h,'Position',[.94 ,pos(2),pos(3)/4,pos(4)])
grid on,title CM_f2-f1,set(gca,'FontSize',7); a=a+1;

% Plot CAPs
h_ax(a)=subplot('position',[0.0446   0.03   0.9429    0.18]); 
% surf([series(:).CAP_ampl],'EdgeAlpha',0)
% set(h_ax(a),'View',[0 -90]),pos=get(h_ax(a),'Position');set(h_ax(a),'Position',[pos(1) ,pos(2),.88,pos(4)]), 
% h = colorbar;pos=get(h,'Position');set(h,'Position',[.94 ,pos(2),pos(3)/4,pos(4)])
% grid on,title CAPs,set(gca,'FontSize',7)

xlabel('no of 6s-presentions'),zoom on
linkaxes(h_ax, 'x')
xlim(h_ax(1),[1 length(series)])
zoom on

