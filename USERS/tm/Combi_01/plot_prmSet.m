function plot_prmSet(results,prmSet,data,number,scrsz)
% TO DO:
% - plot also CM modulation time course of 1st harmonicof f2 and it's
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
if isfield(data,'CAP_unsuppr')
    CAP_dur = length(data(number).CAP_unsuppr);
    figure(get_figure_h([scrsz(1),scrsz(2)+round(scrsz(4)*1/4)-18,round(scrsz(3)/6),...
        round(scrsz(4)/4-40)]));clf
    set(gcf,'Tag','prmSet_fig','NumberTitle','off','Name',['CAP ' results.header.title ' #' num2str(number)])

    plot(linspace(0,1000/f_BT/8,CAP_dur),[data(number).CAP_unsuppr])
    hold on
    plot(linspace(0,2*1000/f_BT,2*8*CAP_dur),[data(number).modCAP; data(number).modCAP])

    yLimit = get(gca,'YLim');
    course = mean(yLimit) + diff(yLimit)/2 * data(number).course_BT_CAP;
    plot(linspace(0,2*1000/f_BT,2*800), [course; course],'k:'),

    % yyaxis right

    plot((5:10:160)*2*1000/f_BT/160,repmat(data(number).CAP_ampl,2,1)/2,'kx')
    plot((5:10:160)*2*1000/f_BT/160,-repmat(data(number).CAP_ampl,2,1)/2,'kx')
    % plot(linspace(0,2*1000/f_BT,2*fs/f_BT),repmat(result(number).CAPmod_course,1,2)/2,'k-')
    % plot(linspace(0,2*1000/f_BT,2*fs/f_BT),-repmat(result(number).CAPmod_course,1,2)/2,'k-')
    CAP_ampl = data(number).CAP_ampl;% + yLimit(1);
    CAPmod_course = spline(5-80:10:75+80,repmat(CAP_ampl,3,1),linspace(0,80,fs/f_BT))';
    plot(linspace(0,2*1000/f_BT,2*fs/f_BT),[CAPmod_course; CAPmod_course]/2,'k-')
    plot(linspace(0,2*1000/f_BT,2*fs/f_BT),-[CAPmod_course; CAPmod_course]/2,'k-')

    xlim([0 2/f_BT*1000]), grid on, zoom yon,
    [~,idx]=max(course(1:end/2));
    set(gca,'XTick',idx*100/3/1600+ [0 .25 .5 .75 1 ]/f_BT*1000,'FontSize',7)
    set(gca,'XTickLabel',{'0' '.25' '.5' '.75' '1' '1.25' '1.5' '1.75' '2'})
    title(['CAP modulation (prmSet: ' num2str(prmSet) ') L_P=' num2str(l2) 'dB L_B=' num2str(l_BT) 'dB, # ' num2str(number)])

    % get delay to display BT phase aligned with that of CAP
    [~,idx_BT]=max(data(number).course_BT(1:fs/f_BT));
    [~,idx_BT_CAP]=max(data(number).course_BT_CAP(1:fs/f_BT));
    circshift_amount = idx_BT_CAP - idx_BT;
else 
    circshift_amount = 0;
end

% plot CM modulation patterns
figure(get_figure_h([scrsz(1),scrsz(2)+round(scrsz(4)*2/4)-18,round(scrsz(3)/6),...
    round(scrsz(4)/4-40)]));clf
set(gcf,'Tag','prmSet_fig','NumberTitle','off','Name',['modCM patterns ' results.header.title ' #' num2str(number)])
plot(x,repmat(circshift(20*log10(abs(data(number).modCM_2f1_f2_course)),circshift_amount),2,1),'b')
hold on;
plot(x,20*log10(abs(data(number).l_unmodCM_2f1_f2)*ones(2/f_BT*fs,1)),'b--')
plot(x,repmat(circshift(20*log10(abs(data(number).modCM_f2_f1_course)),circshift_amount),2,1),'r')
plot(x,20*log10(abs(data(number).l_unmodCM_f2_f1)*ones(2/f_BT*fs,1)),'r--')
% plot(x,repmat(circshift(20*log10(abs(data(number).modCM_2f2_f1_course)),circshift_amount),2,1),'y')
% plot(x,20*log10(abs(data(number).l_unmodCM_2f2_f1)*ones(2/f_BT*fs,1)),'y--')
plot(x,repmat(circshift(20*log10(abs(data(number).modCM_DPf1_course)),circshift_amount),2,1),'Color',[0 .9 0]),
plot(x,20*log10(abs(data(number).l_unmodCM_DPf1)*ones(2/f_BT*fs,1)),'g--','Color',[0 .9 0])
plot(x,repmat(circshift(20*log10(abs(data(number).modCM_DPf2_course)),circshift_amount),2,1),'w')
plot(x,20*log10(abs(data(number).l_unmodCM_DPf2)*ones(2/f_BT*fs,1)),'w--')
plot(x,repmat(circshift(20*log10(abs(data(number).modCM_2f1_course)),circshift_amount),2,1),'m')
plot(x,20*log10(abs(data(number).l_unmodCM_2f1)*ones(2/f_BT*fs,1)),'m--')
plot(x,repmat(circshift(20*log10(abs(data(number).modCM_f2_course)),circshift_amount),2,1),'k')
plot(x,20*log10(abs(data(number).l_unmodCM_f2)*ones(2/f_BT*fs,1)),'k--')

yLimit = get(gca,'YLim');
course = mean(yLimit) + diff(yLimit)/2 * circshift(data(number).course_BT,circshift_amount);
plot(x, [course; course],'k:'),
xlim([0 2/f_BT*1000]), grid on, zoom yon,
[~,idx]=max(course(1:end/2));
set(gca,'XTick',idx*100/3/1600+ [0 .25 .5 .75 1 ]/f_BT*1000,'FontSize',7)
set(gca,'XTickLabel',{'0' '.25' '.5' '.75' '1' '1.25' '1.5' '1.75' '2'})
set(gca,'Color',[.8 .8 .8])
title(['CM modulation (prmSet:' num2str(prmSet) ') L_P=' num2str(l2) 'dB L_B=' num2str(l_BT) 'dB # ' num2str(number)])

% plot OAE modulation patterns
figure(get_figure_h([scrsz(1),scrsz(2)+round(scrsz(4)*3/4)-18,...
    round(scrsz(3)/6),round(scrsz(4)/4-40)]));clf
set(gcf,'Tag','prmSet_fig','NumberTitle','off','Name',['modOAE patterns' results.header.title ' #' num2str(number)])
plot(x,repmat(circshift(20*log10(abs(data(number).mod_2f1_f2_course)),circshift_amount),2,1),'b')
hold on;
plot(x,20*log10(abs(data(number).l_unmod_2f1_f2)*ones(2/f_BT*fs,1)),'b--')
plot(x,repmat(circshift(20*log10(abs(data(number).mod_f2_f1_course)),circshift_amount),2,1),'r')
plot(x,20*log10(abs(data(number).l_unmod_f2_f1)*ones(2/f_BT*fs,1)),'r--')
% plot(x,repmat(circshift(20*log10(abs(data(number).mod_2f2_f1_course)),circshift_amount),2,1),'y')
% plot(x,20*log10(abs(data(number).l_unmod_2f2_f1)*ones(2/f_BT*fs,1)),'y--')
if isfield(data,'modSF_DPf1_course')
    plot(x,repmat(circshift(20*log10(abs(data(number).modSF_DPf1_course)),circshift_amount),2,1),'Color',[0 .9 0])
    plot(x,20*log10(abs(data(number).l_unmod_SF_DPf1)*ones(2/f_BT*fs,1)),'g--','Color',[0 .9 0])
    plot(x,repmat(circshift(20*log10(abs(data(number).modSF_DPf2_course)),circshift_amount),2,1),'w')
    plot(x,20*log10(abs(data(number).l_unmod_SF_DPf2)*ones(2/f_BT*fs,1)),'w--')
    plot(x,repmat(circshift(20*log10(abs(data(number).modSF_f2_course)),circshift_amount),2,1),'k'),
    plot(x,20*log10(abs(data(number).l_unmod_SF_f2)*ones(2/f_BT*fs,1)),'k--')
end
% plot(x,repmat(circshift(20*log10(abs(data(number).mod_2f1_course)),circshift_amount),2,1),'m')
% plot(x,20*log10(abs(data(number).l_unmod_2f1)*ones(2/f_BT*fs,1)),'m--') 
yLimit = get(gca,'YLim');
course = mean(yLimit) + diff(yLimit)/2 * circshift(data(number).course_BT,circshift_amount);
plot(x, [course; course],'k:'),
xlim([0 2/f_BT*1000]), grid on, zoom yon,
[~,idx]=max(course(1:end/2));
set(gca,'XTick',idx*100/3/1600+ [0 .25 .5 .75 1 ]/f_BT*1000,'FontSize',7)
set(gca,'XTickLabel',{'0' '.25' '.5' '.75' '1' '1.25' '1.5' '1.75' '2'})
set(gca,'Color',[.8 .8 .8])
title(['Acoust. modulation (prmSet:' num2str(prmSet) ') L_P=' num2str(l2) 'dB L_B=' num2str(l_BT) 'dB # ' num2str(number)])

% plot(bp_wave), hold on, plot(CAP_samples,zeros(length(CAP_samples),1),'.')

%% plot spectra
if isfield(data,'CAP_unsuppr') % data from ida software (combi.m stimulus)
    H_mod_DPodd = (data(number).H_modDP_1 - data(number).H_modDP_2)/2;
    H_mod_DPeven = (data(number).H_modDP_1 + data(number).H_modDP_2)/2;
    H_mod_DPall = (abs(data(number).H_modDP_1) + abs(data(number).H_modDP_2))/2;
else  % data from OAE.m software
    H_mod_DPodd = data(number).H_modDP_1;
    H_mod_DPeven = data(number).H_modDP_1;
    H_mod_DPall = data(number).H_modDP_1;
end

% plot CM spectrum during modDP
figure(get_figure_h([scrsz(1)+round(scrsz(3)/6),scrsz(2)+round(scrsz(4)*2/4)-18,round(scrsz(3)/6),...
    round(scrsz(4)/4-40)]));clf
set(gcf,'Tag','prmSet_fig','NumberTitle','off','Name',['CM DPspectum ' results.header.title ' #' num2str(number)])
if isfield(data,'CAP_unsuppr') % data from ida software (combi.m stimulus)
    h = stem_60(linspace(0,fs/2-10,fs/20),20*log10(abs(H_mod_DPall(1:fs/20,2))),'.');
    hold on
    if rem(f2,20)~=0  % f1 is multiple of 20
        stem_60((data(number).modlines_2f1_f2-1)*10,20*log10(abs(H_mod_DPodd(data(number).modlines_2f1_f2,2))),'b.')
    else  % This is a quick fix. naming of odd and even reveresed
        stem_60((data(number).modlines_2f1_f2-1)*10,20*log10(abs(H_mod_DPeven(data(number).modlines_2f1_f2,2))),'b.')
    end
    stem_60((data(number).modlines_f2_f1-1)*10,20*log10(abs(H_mod_DPodd(data(number).modlines_f2_f1,2))),'r.')
    stem_60((data(number).modlines_2f2_f1-1)*10,20*log10(abs(data(number).H_modCM_2f2_f1)),'y.')
    stem_60((data(number).modlines_f1-1)*10,20*log10(abs(data(number).H_modCM_DPf1)),'g.')
    stem_60((data(number).modlines_f2-1)*10,20*log10(abs(data(number).H_modCM_DPf2)),'w.')
    stem_60((data(number).modlines_2f1-1)*10,20*log10(abs(H_mod_DPeven(data(number).modlines_2f1,2))),'m.')
else  % data from OAE.m software
    h = stem_60(linspace(0,fs/2-10,fs/10),20*log10(abs(H_mod_DPall(1:fs/10,2))),'.');
    hold on
    stem_60((data(number).modlines_2f1_f2-1)*5,20*log10(abs(H_mod_DPall(data(number).modlines_2f1_f2,2))),'b.')
    stem_60((data(number).modlines_f2_f1-1)*5,20*log10(abs(H_mod_DPodd(data(number).modlines_f2_f1,2))),'r.')
    stem_60((data(number).modlines_2f2_f1-1)*5,20*log10(abs(data(number).H_modCM_2f2_f1)),'y.')
    stem_60((data(number).modlines_f1-1)*5,20*log10(abs(data(number).H_modCM_DPf1)),'g.')
    stem_60((data(number).modlines_f2-1)*5,20*log10(abs(data(number).H_modCM_DPf2)),'w.')
    stem_60((data(number).modlines_2f1-1)*5,20*log10(abs(H_mod_DPeven(data(number).modlines_2f1,2))),'m.')
end
set(h(1),'Color',[.5 .5 .5])
set(h(2),'Color',[.5 .5 .5])
plot(f2,20*log10(abs(data(number).l_unmodCM_DPf2)),'co','MarkerEdgeColor',[.6 .6 .6])
plot(f1,20*log10(abs(data(number).l_unmodCM_DPf1)),'co','MarkerEdgeColor',[.6 1 .6])
plot(2*f1,20*log10(abs(data(number).l_unmodCM_2f1)),'co','MarkerEdgeColor',[1 .6 1])
plot(f2-f1,20*log10(abs(data(number).l_unmodCM_f2_f1)),'co','MarkerEdgeColor',[1 .6 .6])
plot(2*f1-f2,20*log10(abs(data(number).l_unmodCM_2f1_f2)),'co','MarkerEdgeColor',[.6 .6 1])
plot(2*f2-f1,20*log10(abs(data(number).l_unmodCM_2f2_f1)),'co','MarkerEdgeColor',[1 1 .6])
plot([f1 f2],20*log10(abs([data(number).supprCM_f1 data(number).supprCM_f2])),'kx')
grid on,zoom on,axis([1 2.05*f2 -30 60]),set(gca,'FontSize',7),set(gca,'Color',[.8 .8 .8])
title(['CM spectrum modDP (prmSet:' num2str(prmSet) ') L_P=' num2str(l2) 'dB L_B=' num2str(l_BT) 'dB # ' num2str(number)])

% plot CM spectrum during modSF
figure(get_figure_h([scrsz(1)+round(scrsz(3)/6*2),scrsz(2)+round(scrsz(4)*2/4)-18,...
    round(scrsz(3)/6),round(scrsz(4)/4-40)]));clf
set(gcf,'Tag','prmSet_fig','NumberTitle','off','Name',['CM SFspectum ' results.header.title ' #' num2str(number)])
if isfield(data,'CAP_unsuppr') % data from ida software (combi.m stimulus)
    h=stem_60(linspace(0,fs/2-10,fs/20),20*log10(abs(data(number).H_mod_SF(1:fs/20,2))),'.');

    hold on
    stem_60((data(number).modlines_f2-1)*10,20*log10(abs(data(number).H_mod_SF(data(number).modlines_f2,2))),'k.')
else  % data from OAE.m software
    h=stem_60(linspace(0,fs/2-10,fs/10),20*log10(abs(data(number).H_mod_SF(1:fs/10,2))),'.');
    hold on
    stem_60((data(number).modlines_f2-1)*5,20*log10(abs(data(number).H_mod_SF(data(number).modlines_f2,2))),'k.')
end
set(h(1),'Color',[.5 .5 .5])
set(h(2),'Color',[.5 .5 .5])
plot(f2,20*log10(abs(data(number).l_unmodCM_f2)),'co','MarkerEdgeColor',[.6 .6 .6])
grid on,zoom on,axis([1 2.05*f2 -30 60]),set(gca,'FontSize',7),set(gca,'Color',[.8 .8 .8])
title(['CM spectrum modSF (prmSet:' num2str(prmSet) ') L_P=' num2str(l2) 'dB L_B=' num2str(l_BT) 'dB # ' num2str(number)])

% plot acoustic spectrum during modDP
figure(get_figure_h([scrsz(1)+round(scrsz(3)/6),scrsz(2)+round(scrsz(4)*3/4)-18,round(scrsz(3)/6),...
    round(scrsz(4)/4-40)]));clf
set(gcf,'Tag','prmSet_fig','NumberTitle','off','Name',['Acoust DPspectrum ' results.header.title ' #' num2str(number)])
if isfield(data,'CAP_unsuppr') % data from ida software (combi.m stimulus)
    h=stem_60(linspace(0,fs/2-10,fs/20),20*log10(abs(H_mod_DPall(1:fs/20,1))),'.');
    hold on
    if rem(f2,20)~=0  % f1 is multiple of 20
        stem_60((data(number).modlines_2f1_f2-1)*10,20*log10(abs(H_mod_DPodd(data(number).modlines_2f1_f2,1))),'b.')
    else  % This is a quick fix. naming of odd and even reveresed
        stem_60((data(number).modlines_2f1_f2-1)*10,20*log10(abs(H_mod_DPeven(data(number).modlines_2f1_f2,1))),'b.')
    end
    stem_60((data(number).modlines_f2_f1-1)*10,20*log10(abs(H_mod_DPodd(data(number).modlines_f2_f1,1))),'r.')
    stem_60((data(number).modlines_2f2_f1-1)*10,20*log10(abs(data(number).H_mod_2f2_f1)),'y.')
    stem_60((data(number).modlines_f1-1)*10,20*log10(abs(data(number).H_mod_DPf1)),'g.')
    stem_60((data(number).modlines_f2-1)*10,20*log10(abs(data(number).H_mod_DPf2)),'w.')
    stem_60((data(number).modlines_2f1-1)*10,20*log10(abs(H_mod_DPeven(data(number).modlines_2f1,1))),'m.')
else  % data from OAE.m software
    h=stem_60(linspace(0,fs/2-10,fs/10),20*log10(abs(H_mod_DPall(1:fs/10,1))),'.');
    hold on
    stem_60((data(number).modlines_2f1_f2-1)*5,20*log10(abs(H_mod_DPall(data(number).modlines_2f1_f2,1))),'b.')
    stem_60((data(number).modlines_f2_f1-1)*5,20*log10(abs(H_mod_DPodd(data(number).modlines_f2_f1,1))),'r.')
    stem_60((data(number).modlines_2f2_f1-1)*5,20*log10(abs(data(number).H_mod_2f2_f1)),'y.')
    stem_60((data(number).modlines_f1-1)*5,20*log10(abs(data(number).H_mod_DPf1)),'g.')
    stem_60((data(number).modlines_f2-1)*5,20*log10(abs(data(number).H_mod_DPf2)),'w.')
    stem_60((data(number).modlines_2f1-1)*5,20*log10(abs(H_mod_DPeven(data(number).modlines_2f1,1))),'m.')
end
set(h(1),'Color',[.5 .5 .5])
set(h(2),'Color',[.5 .5 .5])
plot([f1 f2 f_BT],[l2+l1 l2 l_BT],'ko')
plot([f1-40 f2-40],[l2+l1+l_suppr l2+l_suppr],'ko')
plot([f1 f2],20*log10(abs([data(number).suppr_f1 data(number).suppr_f2])),'kx')
plot(f2,20*log10(abs(data(number).l_unmod_SF_DPf2)),'co','MarkerEdgeColor',[.6 .6 .6])
plot(f1,20*log10(abs(data(number).l_unmod_SF_DPf1)),'co','MarkerEdgeColor',[.6 1 .6])
plot(f2-f1,20*log10(abs(data(number).l_unmod_f2_f1)),'co','MarkerEdgeColor',[1 .6 .6])
plot(2*f1-f2,20*log10(abs(data(number).l_unmod_2f1_f2)),'co','MarkerEdgeColor',[.6 .6 1])
plot(2*f2-f1,20*log10(abs(data(number).l_unmod_2f2_f1)),'co','MarkerEdgeColor',[1 1 .6])
plot(2*f1,20*log10(abs(data(number).l_unmod_2f1)),'co','MarkerEdgeColor',[1 .6 1]) 
grid on,zoom on,axis([1 2.05*f2 -30 110]),set(gca,'FontSize',7),set(gca,'Color',[.8 .8 .8])
title(['Acoust spectrum modDP (prmSet:' num2str(prmSet) ') L_P=' num2str(l2) 'dB L_B=' num2str(l_BT) 'dB # ' num2str(number)])

% plot acoustic spectrum during modSF
figure(get_figure_h([scrsz(1)+round(scrsz(3)/6*2),scrsz(2)+round(scrsz(4)*3/4)-18,round(scrsz(3)/6),...
    round(scrsz(4)/4-40)]));clf
set(gcf,'Tag','prmSet_fig','NumberTitle','off','Name',['Acoust SFspectum ' results.header.title ' #' num2str(number)])
if isfield(data,'CAP_unsuppr') % data from ida software (combi.m stimulus)
    h=stem_60(linspace(0,fs/2-10,fs/20),20*log10(abs(data(number).H_mod_SF(1:fs/20,1))),'.');
    hold on
    stem_60((data(number).modlines_f2-1)*10,20*log10(abs(data(number).H_mod_f2)),'k.')
else  % data from OAE.m software
    h=stem_60(linspace(0,fs/2-10,fs/10),20*log10(abs(data(number).H_mod_SF(1:fs/10,1))),'.');
    hold on
    stem_60((data(number).modlines_f2-1)*5,20*log10(abs(data(number).H_mod_f2)),'k.')
end
set(h(1),'Color',[.5 .5 .5])
set(h(2),'Color',[.5 .5 .5])
plot(f2,20*log10(abs(data(number).l_unmod_SF_f2)),'co','MarkerEdgeColor',[.6 .6 .6])
grid on,zoom on,axis([1 2.05*f2 -30 110]),set(gca,'FontSize',7),set(gca,'Color',[.8 .8 .8])
title(['Acoust spectrum  modSF (prmSet:' num2str(prmSet) ') L_P=' num2str(l2) 'dB L_B=' num2str(l_BT) 'dB # ' num2str(number)])
