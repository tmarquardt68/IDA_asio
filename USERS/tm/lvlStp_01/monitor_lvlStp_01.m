function result = monitor_combi_01(first,prmSet,results,wave,monitor_settings,plot_flag)
% TO DO:
% - Colour-code phase in time course by marker colour (sub-sampled '.'-markers)
% store mod-lines and compare agerage with avarage modulation time_course
% - Plot time course as "Bode"-diagram.

if ~exist('plot_flag','var')
    plot_flag = 1;
end

latency = 4500;
fs = results.stimulus.sample_rate;
H_mic = results.H_mic;
H_mic = [H_mic(1:end,1) ones(fs/10,1)];
f_BT = results.stimulus.original_parameter_table(prmSet,2);
f1 = results.stimulus.original_parameter_table(prmSet,3);
f2 = results.stimulus.original_parameter_table(prmSet,4);
l1 = results.stimulus.original_parameter_table(prmSet,5);
l2 = results.stimulus.original_parameter_table(prmSet,6);
l_suppr = results.stimulus.original_parameter_table(prmSet,7);
l_BT = results.stimulus.original_parameter_table(prmSet,8);

%% get relevant spectral line numbers

mod_f2_f1(1) = (f2-f1-2*f_BT)/10+1;
mod_f2_f1(2) = (f2-f1-f_BT)/10+1;
mod_f2_f1(3) = (f2-f1)/10+1;
mod_f2_f1(4) = (f2-f1+f_BT)/10+1;
mod_f2_f1(5) = (f2-f1+2*f_BT)/10+1;
result.modlines_f2_f1= mod_f2_f1;

mod_2f1_f2(1) = (2*f1-f2-2*f_BT)/10+1;
mod_2f1_f2(2) = (2*f1-f2-f_BT)/10+1;
mod_2f1_f2(3) = (2*f1-f2)/10+1;
mod_2f1_f2(4) = (2*f1-f2+f_BT)/10+1;
mod_2f1_f2(5) = (2*f1-f2+2*f_BT)/10+1;
result.modlines_2f1_f2 = mod_2f1_f2;

mod_f2(1) = (f2-2*f_BT)/10+1;
mod_f2(2) = (f2-f_BT)/10+1;
mod_f2(3) = (f2)/10+1;
mod_f2(4) = (f2+f_BT)/10+1;
mod_f2(5) = (f2+2*f_BT)/10+1;
result.modlines_f2 = mod_f2;

% get pure speaker primaries in ear canal
H_wave_suppr2 = fft(...
    (wave(16400 + latency+1 : fs/10 + 16400 + latency,:) ...
    - wave(23600 + latency+1 : fs/10 + 23600 + latency,:))/2);
result.l2_suppr = H_wave_suppr2(f2/10+1,1)./H_mic(f2/10+1);
result.supprCM_l2 = H_wave_suppr2(f2/10+1,2);
H_wave_suppr1 = fft(...
    wave(45200 + latency+1 : fs/10 + 45200 + latency,:));
result.l1_suppr = H_wave_suppr1(f1/10+1,1)./H_mic(f1/10+1);
result.supprCM_l1 = H_wave_suppr2(f1/10+1,2);

% get SF0AE level
H_wave_SF = fft(...
    (wave(2000 + latency+1 : fs/10 + 2000 + latency,:) ...
    - wave(9200 + latency+1 : fs/10 + 9200 + latency,:))/2);
result.l_SF = H_wave_SF(f2/10+1,1)./H_mic(f2/10+1)-result.l2_suppr;
result.l_unmodCM_f2 = H_wave_SF(f2/10+1,2);

% get unsuppressed OAE and SF0AE at DP primaries
H_wave_DPodd = fft(...
    (wave(30800 + latency+1 : fs/10 + 30800 + latency,:) ...
    - wave(38000 + latency+1 : fs/10 + 38000 + latency,:))/2);
H_wave_DPeven = fft(...
    (wave(30800 + latency+1 : fs/10 + 30800 + latency,:) ...
    + wave(38000 + latency+1 : fs/10 + 38000 + latency,:))/2);

result.l_2f1_f2 = H_wave_DPodd((2*f1-f2)/10+1,1)./H_mic((2*f1-f2)/10+1);
result.l_unmodCM_2f1_f2 = H_wave_DPodd((2*f1-f2)/10+1,2);
result.l_f2_f1 = H_wave_DPodd((f2-f1)/10+1,1)./H_mic((f2-f1)/10+1);
result.l_unmodCM_f2_f1 = H_wave_DPodd((f2-f1)/10+1,2);
result.l_SF_DPf1 = H_wave_DPeven(f1/10+1,1)./H_mic(f1/10+1)-result.l1_suppr;
result.l_unmodCM_DPf1 = H_wave_DPeven(f1/10+1,2);
result.l_SF_DPf2 = H_wave_DPodd(f2/10+1,1)./H_mic(f2/10+1)-result.l2_suppr;
result.l_unmodCM_DPf2 = H_wave_DPodd(f2/10+1,2);


%% modDP
modDP1_start = 74000 + latency;
avg_wave = zeros(fs/10,2);
for q = 0:7
    avg_wave = avg_wave + wave(q*6*fs/10+modDP1_start+1 : (q*6+1)*fs/10+modDP1_start,:);
end
avg_wave = avg_wave/8;
H_modDP_1=fft(avg_wave)./H_mic;

modDP2_start = 81200 + latency;
avg_wave = zeros(fs/10,2);
for q = 0:7
    avg_wave = avg_wave + wave(q*6*fs/10+modDP2_start+1 : (q*6+1)*fs/10+modDP2_start,:);
end
avg_wave = avg_wave/8;
H_modDP_2=fft(avg_wave)./H_mic;
H_mod_DPodd = (H_modDP_1-H_modDP_2)/2;
H_mod_DPeven = (H_modDP_1+H_modDP_2)/2;
result.l_SF_modDPf2 = H_mod_DPodd(f2/10+1,1) - result.l2_suppr;
result.l_SF_modDPf1 = H_mod_DPeven(f1/10+1,1) - result.l1_suppr;
l = length(H_mod_DPodd);

result.l_BT_ec = H_modDP_2(f_BT/10+1);
H_f_BT = zeros(l,1);
H_f_BT(f_BT/10+1) = result.l_BT_ec;
course = real(ifft(H_f_BT));
result.course_BT = course(1:2/f_BT*fs)/max(course);

%% modSF


modSF1_start = 59600 + latency;
avg_wave = zeros(fs/10,2);
for q = 0:7
    avg_wave = avg_wave + wave(q*6*fs/10+modSF1_start+1 : (q*6+1)*fs/10+modSF1_start,:);
end
avg_wave = avg_wave/8;
H_modSF_1=fft(avg_wave)./H_mic;

modSF2_start = 66800 + latency;
avg_wave = zeros(fs/10,2);
for q = 0:7
    avg_wave = avg_wave + wave(q*6*fs/10+modSF2_start+1 : (q*6+1)*fs/10+modSF2_start,:);
end
avg_wave = avg_wave/8;
H_modSF_2=fft(avg_wave)./H_mic;
H_mod_SF = (H_modSF_1-H_modSF_2)/2;
result.l_modSF = H_mod_SF(f2/10+1) - result.l2_suppr;

% OAE suppression time courses
mod_mask = zeros(l,1);
mod_mask(mod_2f1_f2) = 1;
result.H_mod_2f1_f2 = H_mod_DPodd(mod_2f1_f2(:,1));
course = abs(ifft(H_mod_DPodd(:,1) .* mod_mask));
result.mod_2f1_f2_course = 20*log10(l * course(1:2/f_BT*fs));

mod_mask = zeros(l,1);
mod_mask(mod_f2_f1) = 1;
result.H_mod_f2_f1 = H_mod_DPodd(mod_f2_f1,1);
course = abs(ifft(H_mod_DPodd(:,1) .* mod_mask));
result.mod_f2_f1_course = 20*log10(l * course(1:2/f_BT*fs));

mod_mask = zeros(l,1);
mod_mask(mod_f2) = 1;
H = H_mod_SF(:,1) .* mod_mask;
H(f2/10+1)= result.l_modSF;
result.H_mod_f2 = H(mod_f2);
course = abs(ifft(H));
result.modSF_f2 = 20*log10(l * course(1:2/f_BT*fs));

% CM suppression time courses
mod_mask = zeros(l,1);
mod_mask(mod_2f1_f2) = 1;
result.H_modCM_2f1_f2 = H_mod_DPodd(mod_2f1_f2,2);
course = abs(ifft(H_mod_DPodd(:,2) .* mod_mask));
result.modCM_2f1_f2_course = 20*log10(l * course(1:2/f_BT*fs));

mod_mask = zeros(l,1);
mod_mask(mod_f2_f1) = 1;
result.H_modCM_f2_f1 = H_mod_DPodd(mod_f2_f1,2);
course = abs(ifft(H_mod_DPodd(:,2) .* mod_mask));
result.modCM_f2_f1_course = 20*log10(l * course(1:2/f_BT*fs));

mod_mask = zeros(l,1);
mod_mask(mod_f2) = 1;
result.H_modCM_f2 = H_mod_SF(mod_f2,2);
course = abs(ifft(H_mod_SF(:,2) .* mod_mask));
result.modCM_f2_course = 20*log10(l * course(1:2/f_BT*fs));

% unsuppressed CAP
CAP_latency = latency-300;
CAP_dur = fs/200;
bp_wave = filtfilt(monitor_settings.bp_filter,wave(:,2));
CAP_starts = [2 11 20 29 38 47]*fs/f_BT + CAP_latency;
sign = [1 -1 1 -1 1 -1];
sign = ones(1,6);
avg_wave = zeros(CAP_dur,1); % 20 ms buffer
CAP_samples = [];
for q=1:6
    avg_wave = avg_wave + sign(q)*bp_wave(CAP_starts(q)+1 : CAP_starts(q)+CAP_dur);
    CAP_samples = [CAP_samples CAP_starts(q)+1 : CAP_starts(q)+CAP_dur];
end
result.CAP_unsuppr = avg_wave/6;

% suppressed CAP
modCAP(8*CAP_dur) = 0;
for q=1:8
    modCAP((q-1)*CAP_dur+1:q*CAP_dur) = ...
        (bp_wave((100*(576+(2*q-2)*72+8+q)+1:100*(576+(2*q-2)*72+8+q)+CAP_dur)+CAP_latency) + ...
        bp_wave((100*(576+(2*q-1)*72+8+q)+1:100*(576+(2*q-1)*72+8+q)+CAP_dur)+CAP_latency) +...
        bp_wave((100*(576+(2*q-2+16)*72+8+q)+1:100*(576+(2*q-2+16)*72+8+q)+CAP_dur)+CAP_latency) + ...
        bp_wave((100*(576+(2*q-1+16)*72+8+q)+1:100*(576+(2*q-1+16)*72+8+q)+CAP_dur)+CAP_latency))/4;
    CAP_samples = [CAP_samples [100*(576+(2*q-2)*72+8+q)+1:100*(576+(2*q-2)*72+8+q)+CAP_dur ...
        100*(576+(2*q-1)*72+8+q)+1:100*(576+(2*q-1)*72+8+q)+CAP_dur ...
        100*(576+(2*q-2+16)*72+8+q)+1:100*(576+(2*q-2+16)*72+8+q)+CAP_dur ...
        100*(576+(2*q-1+16)*72+8+q)+1:100*(576+(2*q-1+16)*72+8+q)+CAP_dur]+CAP_latency];
end
result.modCAP = modCAP';

CAP_dur = length(result.CAP_unsuppr);
for q =1:8
    CAP_ampl(q) = max(result.modCAP((q-1)*CAP_dur+1:q*CAP_dur)) - ...
        min(result.modCAP((q-1)*CAP_dur+1:q*CAP_dur));
end
result.CAP_ampl = CAP_ampl;

if ~plot_flag, return, end
%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% plot single measurent
scrsz = get(0,'ScreenSize');

% plot CM spectrum during modDP
figure(get_figure_h([round(scrsz(3)/3.6),round(scrsz(4)*2/4)-18,round(scrsz(3)/3.6),...
    round(scrsz(4)/4-40)]));clf
h = stem_60(linspace(0,fs/2-10,fs/20),20*log10((abs(H_modDP_1(1:fs/20,2))+...
    abs(H_modDP_2(1:fs/20,2)))/2),'.');
set(h(1),'Color',[.5 .5 .5])
set(h(2),'Color',[.5 .5 .5])
hold on
stem_60((mod_2f1_f2-1)*10,20*log10(abs(H_mod_DPodd(mod_2f1_f2,2))),'b.')
stem_60((mod_f2_f1-1)*10,20*log10(abs(H_mod_DPodd(mod_f2_f1,2))),'r.')
plot(f2,20*log10(abs(result.l_unmodCM_DPf2)),'mo')
plot(f1,20*log10(abs(result.l_unmodCM_DPf1)),'mo')
plot(f2,20*log10(abs(result.supprCM_l2)),'go')
plot(f1,20*log10(abs(result.supprCM_l1)),'go')
grid on, zoom on, title('CM spectrum during modDP'), axis([1 24000 -30 60])

% plot acoustic spectrum during modDP
figure(get_figure_h([round(scrsz(3)/3.6),round(scrsz(4)*3/4)-18,round(scrsz(3)/3.6),...
    round(scrsz(4)/4-40)]));clf
h=stem_60(linspace(0,fs/2-10,fs/20),20*log10((abs(H_modDP_1(1:fs/20,1))+...
    abs(H_modDP_2(1:fs/20,1)))/2),'.');
set(h(1),'Color',[.5 .5 .5])
set(h(2),'Color',[.5 .5 .5])
hold on
stem_60((mod_2f1_f2-1)*10,20*log10(abs(H_mod_DPodd(mod_2f1_f2,1))),'b.')
stem_60((mod_f2_f1-1)*10,20*log10(abs(H_mod_DPodd(mod_f2_f1,1))),'r.')
plot([f1 f2 f_BT],[l2+l1 l2 l_BT],'ko')
plot([f1-40 f2-40],[l2+l1+l_suppr l2+l_suppr],'ko')
plot([f1-40 f2-40],[20*log10(abs(H_wave_suppr1(f1/10-3,1)./H_mic(f1/10-3))) ...
    20*log10(abs(H_wave_suppr2(f2/10-3,1)./H_mic(f2/10-3)))],'kx')
plot(f2,20*log10(abs(result.l_SF_DPf2)),'mo')
plot(f1,20*log10(abs(result.l_SF_DPf1)),'mo')
plot(f2,20*log10(abs(result.l_SF_modDPf2)),'go')
plot(f1,20*log10(abs(result.l_SF_modDPf1)),'go')
grid on, zoom on, title('Acoustic spectrum during modDP'), axis([1 24000 -30 110])

% plot CM spectrum during modSF
figure(get_figure_h([round(scrsz(3)/3.6*2),round(scrsz(4)*2/4)-18,...
    round(scrsz(3)/3.6),round(scrsz(4)/4-40)]));clf
h=stem_60(linspace(0,fs/2-10,fs/20),20*log10(abs(H_mod_SF(1:fs/20,2))),'.');
set(h(1),'Color',[.5 .5 .5])
set(h(2),'Color',[.5 .5 .5])
hold on
stem_60((mod_f2-1)*10,20*log10(abs(H_mod_SF(mod_f2,2))),'k.')
plot(f2,20*log10(abs(result.l_unmodCM_f2)),'mo')
grid on, zoom on, title('CM spectrum during modSF'), axis([1 24000 -30 60]),

% plot acoustic spectrum during modSF
figure(get_figure_h([round(scrsz(3)/3.6*2),round(scrsz(4)*3/4)-18,round(scrsz(3)/3.6),...
    round(scrsz(4)/4-40)]));clf
h=stem_60(linspace(0,fs/2-10,fs/20),20*log10(abs(H_mod_SF(1:fs/20,1))),'.');
set(h(1),'Color',[.5 .5 .5])
set(h(2),'Color',[.5 .5 .5])
hold on
stem_60((mod_f2([1 2 4 5])-1)*10,20*log10(abs(H_mod_SF(mod_f2([1 2 4 5]) ,1))),'k.')
plot(f2,20*log10(abs(result.l_SF)),'mo')
plot(f2,20*log10(abs(result.l_modSF)),'go')
grid on, zoom on, title('Acoustic spectrum during modSF'),axis([1 24000 -30 110]), zoom on

% % plot CAPs
% figure(get_figure_h([0,round(scrsz(4)*1/4)-18,round(scrsz(3)/3.6),...
%     round(scrsz(4)/4-40)]));clf
% plot(linspace(0,10,CAP_dur),result.CAP_unsuppr)
% hold on
% plot(linspace(10,90,8*CAP_dur),result.modCAP)
% CAP_ampl2 = 2*(result.CAP_ampl-mean(result.CAP_ampl));
% plot(15:10:85,CAP_ampl2,'kx')
% CAP_mod_curve = spline(5:10:95,[CAP_ampl2(end) CAP_ampl2 CAP_ampl2(1)],linspace(10,90,8*CAP_dur));
% plot(linspace(10,90,8*CAP_dur),CAP_mod_curve,'k-')
% grid on, title(['CAP - L_B_T= ' num2str(l_BT) 'dB, L2=' num2str(l2) 'dB @ ' num2str(f2) 'Hz'])

plot_prmSet(results, prmSet, result)

%% plot time series of relevant spectral line levels
h_series = get_figure_h([1385 1 540 850]);
figure(h_series),
set(h_series,'Name','Series')
if first==1
    clf
    modDPlinesSeries = [];
    modSFlinesSeries = [];
    setappdata(h_series,'modDPlinesSeries',modDPlinesSeries)
    setappdata(h_series,'modSFlinesSeries',modSFlinesSeries)
end
modDPlinesSeries = getappdata(h_series,'modDPlinesSeries');
modSFlinesSeries = getappdata(h_series,'modSFlinesSeries');

modDPlinesSeries = [modDPlinesSeries [H_mod_DPodd([mod_2f1_f2 mod_f2_f1],1);...
    H_mod_DPodd([mod_2f1_f2 mod_f2_f1],2)]];
modSFlinesSeries = [modSFlinesSeries [H_mod_SF(mod_f2,1);result.l_SF;H_mod_SF(mod_f2,2)]];
modSFlinesSeries(3,end) = result.l_modSF;

setappdata(h_series,'modDPlinesSeries',modDPlinesSeries)
setappdata(h_series,'modSFlinesSeries',modSFlinesSeries)
plot_time_series(modDPlinesSeries,modSFlinesSeries)








