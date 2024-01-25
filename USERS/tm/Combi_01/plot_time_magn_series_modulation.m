function plot_time_magn_series_modulation(results,series)
n = length(series);
fs = results.stimulus.sample_rate;
interval_order=results.stimulus.interval_order';
f_BT(2) = results.stimulus.original_parameter_table(interval_order(2),2);
if isfield(results.data,'CAP_unsuppr')
    f_BT = zeros(1,n);
    circshift_amount = f_BT;
for q=2 % 1:n Loop on case f_BT is chaning in this series
    f_BT(q) = results.stimulus.original_parameter_table(interval_order(q),2);
    [~,idx_BT]=max(series(q).course_BT(1:fs/f_BT(q)));
    [~,idx_BT_CAP]=max(series(q).course_BT_CAP(1:fs/f_BT(q)));
    circshift_amount(q) = idx_BT_CAP - idx_BT;
    [~,idx]=max(series(q).course_BT(1:fs/f_BT(q)));
    circshift_amount(q) = -circshift_amount(q) + idx;
    circshift_amount_CAP(q) = -idx_BT_CAP;
end
else 
    [~,idx_BT]=max(series(2).course_BT(1:fs/f_BT(2)));
    circshift_amount(2) = -idx_BT;
end
a=1;
% mod_2f1_f2_course
h_ax(a)=subplot('position',[0.0446    0.88    0.9429    0.1]);
course = circshift(20*log10(abs([series(:).mod_2f1_f2_course])),circshift_amount(2));
surf([course course(:,end)],'EdgeAlpha',0)
pos=get(h_ax(a),'Position');set(h_ax(a),'Position',[pos(1) ,pos(2),.88,pos(4)]), 
set(h_ax(a),'View',[0 -90],'YColor',[0 0 1],'YTick',[0 .25 .5 .75 1 ]*fs/f_BT(2))
set(gca,'YTickLabel',{'0' '.25' '.5' '.75' '1'})
h = colorbar;pos=get(h,'Position');set(h,'Position',[.94 ,pos(2),pos(3)/4,pos(4)])
grid on,title('2f1-f1'),set(gca,'FontSize',7); a=a+1;

% modCM_2f1_f2_course
h_ax(a)=subplot('position',[0.0446    0.76    0.9429    0.1]);
course = circshift(20*log10(abs([series(:).modCM_2f1_f2_course])),circshift_amount(2));
surf([course course(:,end)],'EdgeAlpha',0)
pos=get(h_ax(a),'Position');set(h_ax(a),'Position',[pos(1) ,pos(2),.88,pos(4)]), 
set(h_ax(a),'View',[0 -90],'YColor',[0 0 1],'YTick',[0 .25 .5 .75 1 ]*fs/f_BT(2))
set(gca,'YTickLabel',{'0' '.25' '.5' '.75' '1'})
h = colorbar;pos=get(h,'Position');set(h,'Position',[.94 ,pos(2),pos(3)/4,pos(4)])
grid on,title('CM 2f1-f2'),set(gca,'FontSize',7); a=a+1;

% mod_f2_f1_course
h_ax(a)=subplot('position',[0.0446    0.64    0.9429    0.1]);
course = circshift(20*log10(abs([series(:).mod_f2_f1_course])),circshift_amount(2));
surf([course course(:,end)],'EdgeAlpha',0)
pos=get(h_ax(a),'Position');set(h_ax(a),'Position',[pos(1) ,pos(2),.88,pos(4)]), 
set(h_ax(a),'View',[0 -90],'YColor',[1 0 0],'YTick',[0 .25 .5 .75 1 ]*fs/f_BT(2))
set(gca,'YTickLabel',{'0' '.25' '.5' '.75' '1'})
h = colorbar;pos=get(h,'Position');set(h,'Position',[.94 ,pos(2),pos(3)/4,pos(4)])
grid on,title('f2-f1'),set(gca,'FontSize',7); a=a+1;

% modCM_f2_f1_course
h_ax(a)=subplot('position',[0.0446   0.52   0.9429    0.1]);
course = circshift(20*log10(abs([series(:).modCM_f2_f1_course])),circshift_amount(2));
surf([course course(:,end)],'EdgeAlpha',0)
pos=get(h_ax(a),'Position');set(h_ax(a),'Position',[pos(1) ,pos(2),.88,pos(4)]), 
set(h_ax(a),'View',[0 -90],'YColor',[1 0 0],'YTick',[0 .25 .5 .75 1 ]*fs/f_BT(2))
set(gca,'YTickLabel',{'0' '.25' '.5' '.75' '1'})
h = colorbar;pos=get(h,'Position');set(h,'Position',[.94 ,pos(2),pos(3)/4,pos(4)])
grid on,title('CM f2-f1'),set(gca,'FontSize',7); a=a+1;

% modSF_f2_course 
h_ax(a)=subplot('position',[0.0446    0.40    0.9429    0.1]);
course = circshift(20*log10(abs([series(:).modSF_f2_course])),circshift_amount(2));
surf([course course(:,end)],'EdgeAlpha',0)
pos=get(h_ax(a),'Position');set(h_ax(a),'Position',[pos(1) ,pos(2),.88,pos(4)]), 
set(h_ax(a),'View',[0 -90],'YColor',[0 0 0],'YTick',[0 .25 .5 .75 1 ]*fs/f_BT(2))
set(gca,'YTickLabel',{'0' '.25' '.5' '.75' '1'})
h = colorbar;pos=get(h,'Position');set(h,'Position',[.94 ,pos(2),pos(3)/4,pos(4)])
grid on,title('f2'),set(gca,'FontSize',7); a=a+1;

% modCM_f2_course level
h_ax(a)=subplot('position',[0.0446   0.28   0.9429    0.1]);
course = circshift(20*log10(abs([series(:).modCM_f2_course])),circshift_amount(2));
surf([course course(:,end)],'EdgeAlpha',0)
pos=get(h_ax(a),'Position');set(h_ax(a),'Position',[pos(1) ,pos(2),.88,pos(4)]), 
set(h_ax(a),'View',[0 -90],'YColor',[0 0 0],'YTick',[0 .25 .5 .75 1 ]*fs/f_BT(2))
set(gca,'YTickLabel',{'0' '.25' '.5' '.75' '1'})
h = colorbar;pos=get(h,'Position');set(h,'Position',[.94 ,pos(2),pos(3)/4,pos(4)])
grid on,title('CM f2'),set(gca,'FontSize',7); a=a+1;

% modCM_2f1_course level
h_ax(a)=subplot('position',[0.0446   0.16   0.9429    0.1]);
course = circshift(20*log10(abs([series(:).modCM_2f1_course])),circshift_amount(2));
surf([course course(:,end)],'EdgeAlpha',0)
pos=get(h_ax(a),'Position');set(h_ax(a),'Position',[pos(1) ,pos(2),.88,pos(4)]), 
set(h_ax(a),'View',[0 -90],'YColor',[1 0 1],'YTick',[0 .25 .5 .75 1 ]*fs/f_BT(2))
set(gca,'YTickLabel',{'0' '.25' '.5' '.75' '1'})
h = colorbar;pos=get(h,'Position');set(h,'Position',[.94 ,pos(2),pos(3)/4,pos(4)])
grid on,title('CM 2f1'),set(gca,'FontSize',7); a=a+1;

% CAP_ampl
if isfield(results.data,'CAP_unsuppr')
    h_ax(a)=subplot('position',[0.0446   0.03   0.9429    0.1]);
    for q=1:n
        CAP_ampl = series(q).CAP_ampl;
        CAPmod_course(:,q) = spline(5-80:10:75+80,[CAP_ampl; CAP_ampl; CAP_ampl],linspace(0,80,fs/f_BT(2)))';
    end
    course = circshift(CAPmod_course,circshift_amount_CAP(2))*1000;
    % course = circshift([series(:).CAPmod_course],circshift_amount(2))*1000;
    pos=get(h_ax(a),'Position');set(h_ax(a),'Position',[pos(1) ,pos(2),.88,pos(4)]),
    if isfield(series,'CAP_ampl')
        surf([course course(:,end)],'EdgeAlpha',0)
        h = colorbar;pos=get(h,'Position');
        set(h,'Position',[.94 ,pos(2),pos(3)/4,pos(4)])
    end
    set(h_ax(a),'View',[0 -90],'YColor',[0 0 1],'YTick',[0 .25 .5 .75 1 ]*fs/f_BT(2))
    set(gca,'YTickLabel',{'0' '.25' '.5' '.75' '1'})
    grid on,title('CAPs'), set(gca,'FontSize',7,'Tag','axes_CAP')
end

xlabel('no of 6s-presentions')
linkaxes(h_ax, 'x')
xlim(h_ax(1),[1 length(series)+1])

