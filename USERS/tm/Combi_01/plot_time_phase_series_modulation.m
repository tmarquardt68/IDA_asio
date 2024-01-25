function plot_time_phase_series_modulation(results,series)
n = length(series);
fs = results.stimulus.sample_rate;
interval_order=results.stimulus.interval_order';
f_BT = zeros(1,n);
circshift_amount = f_BT;
for q=2 % 1:n  % !!! loop will be required if f_BT is not constangt
    f_BT(q) = results.stimulus.original_parameter_table(interval_order(q),2);
    [~,idx_BT]=max(series(q).course_BT(1:fs/f_BT(q)));
    [~,idx_BT_CAP]=max(series(q).course_BT_CAP(1:fs/f_BT(q)));
    circshift_amount(q) = idx_BT_CAP - idx_BT;
    [~,idx]=max(series(q).course_BT(1:fs/f_BT(q)));
    circshift_amount(q) = -circshift_amount(q) + idx;
end

a=1;
% mod_2f1_f2_course
h_ax(a)=subplot('position',[0.0446    0.88    0.9429    0.1]);
surf(circshift(angle([series(:).mod_2f1_f2_course]/series(1).l_unmod_2f1_f2),circshift_amount(2))/2/pi,'EdgeAlpha',0)
pos=get(h_ax(a),'Position');set(h_ax(a),'Position',[pos(1) ,pos(2),.88,pos(4)]), 
set(h_ax(a),'View',[0 -90],'YColor',[0 0 1],'YTick',[0 .25 .5 .75 1 ]*fs/f_BT(2))
set(gca,'YTickLabel',{'0' '.25' '.5' '.75' '1'})
h = colorbar;pos=get(h,'Position');set(h,'Position',[.94 ,pos(2),pos(3)/4,pos(4)])
grid on,title('2f1-f1'),set(gca,'FontSize',7); a=a+1;

% modCM_2f1_f2_course
h_ax(a)=subplot('position',[0.0446    0.76    0.9429    0.1]);
surf(circshift(angle([series(:).modCM_2f1_f2_course]/series(1).l_unmodCM_2f2_f1),circshift_amount(2))/2/pi,'EdgeAlpha',0)
pos=get(h_ax(a),'Position');set(h_ax(a),'Position',[pos(1) ,pos(2),.88,pos(4)]), 
set(h_ax(a),'View',[0 -90],'YColor',[0 0 1],'YTick',[0 .25 .5 .75 1 ]*fs/f_BT(2))
set(gca,'YTickLabel',{'0' '.25' '.5' '.75' '1'})
h = colorbar;pos=get(h,'Position');set(h,'Position',[.94 ,pos(2),pos(3)/4,pos(4)])
grid on,title('CM 2f1-f2'),set(gca,'FontSize',7); a=a+1;

% mod_f2_f1_course
h_ax(a)=subplot('position',[0.0446    0.64    0.9429    0.1]);
surf(circshift(angle([series(:).mod_f2_f1_course]/series(1).l_unmod_f2_f1),circshift_amount(2))/2/pi,'EdgeAlpha',0)
pos=get(h_ax(a),'Position');set(h_ax(a),'Position',[pos(1) ,pos(2),.88,pos(4)]), 
set(h_ax(a),'View',[0 -90],'YColor',[1 0 0],'YTick',[0 .25 .5 .75 1 ]*fs/f_BT(2))
set(gca,'YTickLabel',{'0' '.25' '.5' '.75' '1'})
h = colorbar;pos=get(h,'Position');set(h,'Position',[.94 ,pos(2),pos(3)/4,pos(4)])
grid on,title('f2-f1'),set(gca,'FontSize',7); a=a+1;

% modCM_f2_f1_course
h_ax(a)=subplot('position',[0.0446   0.52   0.9429    0.1]);
surf(circshift(angle([series(:).modCM_f2_f1_course]/series(1).l_unmodCM_f2_f1),circshift_amount(2))/2/pi,'EdgeAlpha',0)
pos=get(h_ax(a),'Position');set(h_ax(a),'Position',[pos(1) ,pos(2),.88,pos(4)]), 
set(h_ax(a),'View',[0 -90],'YColor',[1 0 0],'YTick',[0 .25 .5 .75 1 ]*fs/f_BT(2))
set(gca,'YTickLabel',{'0' '.25' '.5' '.75' '1'})
h = colorbar;pos=get(h,'Position');set(h,'Position',[.94 ,pos(2),pos(3)/4,pos(4)])
grid on,title('CM f2-f1'),set(gca,'FontSize',7); a=a+1;

% modSF_f2_course 
h_ax(a)=subplot('position',[0.0446    0.40    0.9429    0.1]);
surf(circshift(angle([series(:).modSF_f2_course]/series(1).l_unmod_SF_f2),circshift_amount(2))/2/pi,'EdgeAlpha',0)
pos=get(h_ax(a),'Position');set(h_ax(a),'Position',[pos(1) ,pos(2),.88,pos(4)]), 
set(h_ax(a),'View',[0 -90],'YColor',[0 0 0],'YTick',[0 .25 .5 .75 1 ]*fs/f_BT(2))
set(gca,'YTickLabel',{'0' '.25' '.5' '.75' '1'})
h = colorbar;pos=get(h,'Position');set(h,'Position',[.94 ,pos(2),pos(3)/4,pos(4)])
grid on,title('f2'),set(gca,'FontSize',7); a=a+1;

% modCM_f2_course level
h_ax(a)=subplot('position',[0.0446   0.28   0.9429    0.1]);
surf(circshift(angle([series(:).modCM_f2_course]/series(1).l_unmodCM_f2),circshift_amount(2))/2/pi,'EdgeAlpha',0)
pos=get(h_ax(a),'Position');set(h_ax(a),'Position',[pos(1) ,pos(2),.88,pos(4)]), 
set(h_ax(a),'View',[0 -90],'YColor',[0 0 0],'YTick',[0 .25 .5 .75 1 ]*fs/f_BT(2))
set(gca,'YTickLabel',{'0' '.25' '.5' '.75' '1'})
h = colorbar;pos=get(h,'Position');set(h,'Position',[.94 ,pos(2),pos(3)/4,pos(4)])
grid on,title('CM f2'),set(gca,'FontSize',7); a=a+1;

% CAP_ampl
h_ax(a)=subplot('position',[0.0446   0.03   0.9429    0.1]); 
pos=get(h_ax(a),'Position');set(h_ax(a),'Position',[pos(1) ,pos(2),.88,pos(4)]),
if isfield(series,'CAP_ampl')
    surf(circshift([series(:).CAP_ampl],circshift_amount(2)),'EdgeAlpha',0)
    h = colorbar;pos=get(h,'Position');
    set(h,'Position',[.94 ,pos(2),pos(3)/4,pos(4)])
end
set(h_ax(a),'View',[0 -90],'YColor',[0 0 1],'YTick',[0 .25 .5 .75 1 ]*fs/f_BT(2))
set(gca,'YTickLabel',{'0' '.25' '.5' '.75' '1'})
grid on,title('CAPs'), set(gca,'FontSize',7,'Tag','axes_CAP')

xlabel('no of 6s-presentions')
linkaxes(h_ax, 'x')
xlim(h_ax(1),[1 length(series)])

