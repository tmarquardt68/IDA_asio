function result = monitor_combi_01(first,prmSet,results,wave,monitor_settings,plot_flag)
% TO DO:
% - "l_unmod*" => "unmod*" , "l_mod*" to "mod*" , "l_SF_modDPf1" =>
% "modSF_DPf1" 

monitor = 1; % chose monitor
scrsz = get(0,'MonitorPositions');

if ~exist('plot_flag','var')
    plot_flag = 1;
end

latency = 4500;
fs = results.stimulus.sample_rate;
H_mic = results.H_mic;
% H_mic = [H_mic(1:fs/10,1) ones(fs/10,1)];  !!! uncomment with proper
% calibration (first when started with GP experiments in April 2023)
% H_mic = ones(fs/10,2);
if length(H_mic)==fs/10
    H_mic = [H_mic ones(fs/10,1)];
else
    H_mic = [H_mic(1:2:end,1) ones(fs/10,1)];
end
f_BT = results.stimulus.original_parameter_table(prmSet,2);
f1 = results.stimulus.original_parameter_table(prmSet,3);
f2 = results.stimulus.original_parameter_table(prmSet,4);
l = fs/10;

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

mod_2f2_f1(1) = (2*f2-f1-2*f_BT)/10+1;
mod_2f2_f1(2) = (2*f2-f1-f_BT)/10+1;
mod_2f2_f1(3) = (2*f2-f1)/10+1;
mod_2f2_f1(4) = (2*f2-f1+f_BT)/10+1;
mod_2f2_f1(5) = (2*f2-f1+2*f_BT)/10+1;
result.modlines_2f2_f1 = mod_2f2_f1;

mod_f2(1) = (f2-2*f_BT)/10+1;
mod_f2(2) = (f2-f_BT)/10+1;
mod_f2(3) = (f2)/10+1;
mod_f2(4) = (f2+f_BT)/10+1;
mod_f2(5) = (f2+2*f_BT)/10+1;
result.modlines_f2 = mod_f2;

mod_f1(1) = (f1-2*f_BT)/10+1;
mod_f1(2) = (f1-f_BT)/10+1;
mod_f1(3) = (f1)/10+1;
mod_f1(4) = (f1+f_BT)/10+1;
mod_f1(5) = (f1+2*f_BT)/10+1;
result.modlines_f1 = mod_f1;

mod_2f1(1) = (2*f1-2*f_BT)/10+1;
mod_2f1(2) = (2*f1-f_BT)/10+1;
mod_2f1(3) = (2*f1)/10+1;
mod_2f1(4) = (2*f1+f_BT)/10+1;
mod_2f1(5) = (2*f1+2*f_BT)/10+1;
result.modlines_2f1= mod_2f1;

%% CAP
% unsuppressed CAP
CAPwindowStart_latency = latency-300;
CAP_dur = fs/200;
bp_wave = filtfilt(monitor_settings.bp_filter,wave(:,2));
CAP_starts = [2 11 20 29 38 47]*fs/f_BT + CAPwindowStart_latency;
for q=6:-1:1
    CAP_unsuppr(q,:) = bp_wave(CAP_starts(q)+1 : CAP_starts(q)+CAP_dur);
end
result.CAP_unsuppr = mean(CAP_unsuppr);

% suppressed CAP
for q=8:-1:1
    modCAP((q-1)*CAP_dur+1:q*CAP_dur) = ...
        (bp_wave((100*(576+(2*q-2)*72+8+q)+1:100*(576+(2*q-2)*72+8+q)+CAP_dur)+CAPwindowStart_latency) + ...
        bp_wave((100*(576+(2*q-1)*72+8+q)+1:100*(576+(2*q-1)*72+8+q)+CAP_dur)+CAPwindowStart_latency) +...
        bp_wave((100*(576+(2*q-2+16)*72+8+q)+1:100*(576+(2*q-2+16)*72+8+q)+CAP_dur)+CAPwindowStart_latency) + ...
        bp_wave((100*(576+(2*q-1+16)*72+8+q)+1:100*(576+(2*q-1+16)*72+8+q)+CAP_dur)+CAPwindowStart_latency))/4;
end
result.modCAP = modCAP';

for q =8:-1:1
    CAP_ampl(q) = max(result.modCAP((q-1)*CAP_dur+1:q*CAP_dur)) - ...
        min(result.modCAP((q-1)*CAP_dur+1:q*CAP_dur));
end
result.CAP_ampl = CAP_ampl';
result.CAPmod_course = spline(5-80:10:75+80,repmat(CAP_ampl,1,3),linspace(0,80,fs/f_BT))';

% get time course of BT related to modCAP
for q=8:-1:1
    BT_ampl(:,q) = ...
        [wave((100*(576+(2*q-2)*72+8+q)+1)+CAPwindowStart_latency,1) ; ...
        wave((100*(576+(2*q-1)*72+8+q)+1)+CAPwindowStart_latency,1) ;...
        wave((100*(576+(2*q-2+16)*72+8+q)+1)+CAPwindowStart_latency,1) ; ...
        wave((100*(576+(2*q-1+16)*72+8+q)+1)+CAPwindowStart_latency,1)];
end
H_f_BT = zeros(l,1);
tmp = fft(mean(BT_ampl));
% delay so that time when BT amplitude is read
% (CAPwindowStart_latency) is at f2 onset and also aligned with CAP peak
H_f_BT(f_BT/10+1) = tmp(2)*exp(-1i*0.125e-3*f_BT*2*pi);
course = real(ifft(H_f_BT));
result.course_BT_CAP = course(1:fs/f_BT)/max(course);

%% biasing tone OFF %%
% get pure speaker primaries in ear canal during SFOAE suppression
if rem(f2,20)~=0  % f1 is multiple of 20
    avg_wave = (wave(16400 + latency+1 : fs/10 + 16400 + latency,:) ...
        - wave(23600 + latency+1 : fs/10 + 23600 + latency,:))/2;
    H_wave_suppr2 = fft(avg_wave)./H_mic;
    result.suppr_f2 = H_wave_suppr2(f2/10+1,1);
    result.supprCM_f2 = H_wave_suppr2(f2/10+1,2);

    avg_wave =  wave(45200 + latency+1 : fs/10 + 45200 + latency,:);
    H_wave_suppr1 = fft(avg_wave)./H_mic;
    result.suppr_f1 = H_wave_suppr1(f1/10+1,1);
    result.supprCM_f1 = H_wave_suppr1(f1/10+1,2);
    result.supprCM_2f1 = H_wave_suppr1(f1/10+1,2);

    % get unmodulated SF0AE level
    avg_wave = (wave(2000 + latency+1 : fs/10 + 2000 + latency,:) ...
        - wave(9200 + latency+1 : fs/10 + 9200 + latency,:))/2;
    H_wave_SF = fft(avg_wave)./H_mic;
    result.l_unmod_SF_f2 = H_wave_SF(f2/10+1,1)-result.suppr_f2;
    result.l_unmodCM_f2 = H_wave_SF(f2/10+1,2);
    result.l_unmod_f2 = H_wave_SF(f2/10+1,1);

    % get unmodulated DPOAE and SF0AE at DP primaries
    avg_wave = (wave(30800 + latency+1 : fs/10 + 30800 + latency,:) ...
        - wave(38000 + latency+1 : fs/10 + 38000 + latency,:))/2;
    H_wave_unmod_DPodd = fft(avg_wave)./H_mic;

    avg_wave = (wave(30800 + latency+1 : fs/10 + 30800 + latency,:) ...
        + wave(38000 + latency+1 : fs/10 + 38000 + latency,:))/2;
    H_wave_unmod_DPeven = fft(avg_wave)./H_mic;

    result.l_unmod_2f1_f2 = H_wave_unmod_DPodd((2*f1-f2)/10+1,1);
    result.l_unmodCM_2f1_f2 = H_wave_unmod_DPodd((2*f1-f2)/10+1,2);
    result.l_unmod_f2_f1 = H_wave_unmod_DPodd((f2-f1)/10+1,1);
    result.l_unmodCM_f2_f1 = H_wave_unmod_DPodd((f2-f1)/10+1,2);
    result.l_unmod_2f2_f1 = H_wave_unmod_DPeven((2*f2-f1)/10+1,1);
    result.l_unmodCM_2f2_f1 = H_wave_unmod_DPeven((2*f2-f1)/10+1,2);
    result.l_unmod_SF_DPf1 = H_wave_unmod_DPeven(f1/10+1,1)-result.suppr_f1;
    result.l_unmodCM_DPf1 = H_wave_unmod_DPeven(f1/10+1,2);
    result.l_unmod_SF_DPf2 = H_wave_unmod_DPodd(f2/10+1,1)-result.suppr_f2;
    result.l_unmodCM_DPf2 = H_wave_unmod_DPodd(f2/10+1,2);
    result.l_unmod_2f1 = H_wave_unmod_DPeven(2*f1/10+1,1);
    result.l_unmodCM_2f1 = H_wave_unmod_DPeven(2*f1/10+1,2);

    %% biasing tone ON %%
    %% simulataneously presented primaries (modDP)
    modDP1_start = 74000 + latency;
    avg_wave = zeros(fs/10,2);
    for q = 0:7
        avg_wave = avg_wave + wave(q*6*fs/10+modDP1_start+1 : (q*6+1)*fs/10+modDP1_start,:);
    end
    avg_wave = avg_wave/8;
    result.H_modDP_1=fft(avg_wave)./H_mic;

    modDP2_start = 81200 + latency;
    avg_wave = zeros(fs/10,2);
    for q = 0:7
        avg_wave = avg_wave + wave(q*6*fs/10+modDP2_start+1 : (q*6+1)*fs/10+modDP2_start,:);
    end
    avg_wave = avg_wave/8;
    result.H_modDP_2=fft(avg_wave)./H_mic;

    H_mod_DPodd = (result.H_modDP_1-result.H_modDP_2)/2;
    H_mod_DPeven = (result.H_modDP_1+result.H_modDP_2)/2;

    % BT time course
    result.l_BT_ec = result.H_modDP_2(f_BT/10+1);
    H_f_BT = zeros(l,1);
    H_f_BT(f_BT/10+1) = result.l_BT_ec;
    course = real(ifft(H_f_BT));
    result.course_BT = course(1:fs/f_BT)/max(course);

    % 2f1-f2 suppression time courses
    mod_mask = zeros(l,1);
    mod_mask(mod_2f1_f2) = 1;

    result.H_mod_2f1_f2 = H_mod_DPodd(mod_2f1_f2,1);
    course = l * ifft(H_mod_DPodd(:,1) .* mod_mask);
    result.mod_2f1_f2_course = course(1:fs/f_BT);

    result.H_modCM_2f1_f2 = H_mod_DPodd(mod_2f1_f2,2);
    course = l * ifft(H_mod_DPodd(:,2) .* mod_mask);
    result.modCM_2f1_f2_course = course(1:fs/f_BT);

    % f2-f1 suppression time courses
    mod_mask = zeros(l,1);
    mod_mask(mod_f2_f1) = 1;

    result.H_mod_f2_f1 = H_mod_DPodd(mod_f2_f1,1);
    course = l * ifft(H_mod_DPodd(:,1) .* mod_mask);
    result.mod_f2_f1_course = course(1:fs/f_BT);

    result.H_modCM_f2_f1 = H_mod_DPodd(mod_f2_f1,2);
    course = l * ifft(H_mod_DPodd(:,2) .* mod_mask);
    result.modCM_f2_f1_course = course(1:fs/f_BT);

    % 2f2-f1 suppression time courses
    mod_mask = zeros(l,1);
    mod_mask(mod_2f2_f1) = 1;

    result.H_modCM_2f2_f1 = H_mod_DPeven(mod_2f2_f1,2);
    course = l * ifft(H_mod_DPeven(:,2) .* mod_mask);
    result.modCM_2f2_f1_course = course(1:fs/f_BT);

    result.H_mod_2f2_f1 = H_mod_DPeven(mod_2f2_f1,1);
    course = l * ifft(H_mod_DPeven(:,1) .* mod_mask);
    result.mod_2f2_f1_course = course(1:fs/f_BT);

    % 2f1 suppression time courses
    mod_mask = zeros(l,1);
    mod_mask(mod_2f1) = 1;

    result.H_modCM_2f1 = H_mod_DPeven(mod_2f1,2);
    course = l * ifft(H_mod_DPeven(:,2) .* mod_mask);
    result.modCM_2f1_course = course(1:fs/f_BT);

    result.H_mod_2f1 = H_mod_DPeven(mod_2f1,1);
    course = l * ifft(H_mod_DPeven(:,1) .* mod_mask);
    result.mod_2f1_course = course(1:fs/f_BT);

    % time courses of f1-SFOAE produced during simulataneously presented primaries and biasing tone
    result.l_SF_modDPf1 = H_mod_DPeven(f1/10+1,1) - result.suppr_f1;
    mod_mask = zeros(l,1);
    mod_mask(mod_f1) = 1;

    result.H_modCM_DPf1 = H_mod_DPeven(mod_f1,2);
    course = l * ifft(H_mod_DPeven(:,2) .* mod_mask);
    result.modCM_DPf1_course = course(1:fs/f_BT);

    H = H_mod_DPeven(:,1) .* mod_mask;
    H(f1/10+1)= result.l_SF_modDPf1;
    result.H_mod_DPf1 = H(mod_f1);
    course = l * ifft(H);
    result.modSF_DPf1_course = course(1:fs/f_BT);

    % time courses of f2-SFOAE produced during simulataneously presented primaries and biasing tone
    result.l_SF_modDPf2 = H_mod_DPodd(f2/10+1,1) - result.suppr_f2;
    mod_mask = zeros(l,1);
    mod_mask(mod_f2) = 1;

    result.H_modCM_DPf2 = H_mod_DPodd(mod_f2,2);
    course = l * ifft(H_mod_DPodd(:,2) .* mod_mask);
    result.modCM_DPf2_course = course(1:fs/f_BT);

    H = H_mod_DPodd(:,1) .* mod_mask;
    H(f2/10+1)= result.l_SF_modDPf2;
    result.H_mod_DPf2 = H(mod_f2);
    course = l * ifft(H);
    result.modSF_DPf2_course = course(1:fs/f_BT);


    %% only f2 primary presented (modSF)
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
    result.H_mod_SF = (H_modSF_1-H_modSF_2)/2;
    result.l_modSF_f2 = result.H_mod_SF(f2/10+1) - result.suppr_f2;

    % SFOAE time course
    mod_mask = zeros(l,1);
    mod_mask(mod_f2) = 1;

    result.H_modCM_f2 = result.H_mod_SF(mod_f2,2);
    course = l * ifft(result.H_mod_SF(:,2) .* mod_mask);
    result.modCM_f2_course = course(1:fs/f_BT);

    H = result.H_mod_SF(:,1) .* mod_mask;
    H(f2/10+1)= result.l_modSF_f2;
    result.H_mod_f2 = H(mod_f2);
    course = l * ifft(H);
    result.modSF_f2_course = course(1:fs/f_BT);

else % f2 is multiple of 20
    warning('TM warns: f2 is multiple of 20. CAP might not be sum of in phase and inverted phase onset!')
    avg_wave = (wave(16400 + latency+1 : fs/10 + 16400 + latency,:) ...
        + wave(23600 + latency+1 : fs/10 + 23600 + latency,:))/2;
    H_wave_suppr2 = fft(avg_wave)./H_mic;
    result.suppr_f2 = H_wave_suppr2(f2/10+1,1);
    result.supprCM_f2 = H_wave_suppr2(f2/10+1,2);

    avg_wave =  wave(45200 + latency+1 : fs/10 + 45200 + latency,:);
    H_wave_suppr1 = fft(avg_wave)./H_mic;
    result.suppr_f1 = H_wave_suppr1(f1/10+1,1);
    result.supprCM_f1 = H_wave_suppr1(f1/10+1,2);
    result.supprCM_2f1 = H_wave_suppr1(f1/10+1,2);

    % get unmodulated SF0AE level
    avg_wave = (wave(2000 + latency+1 : fs/10 + 2000 + latency,:) ...
        + wave(9200 + latency+1 : fs/10 + 9200 + latency,:))/2;
    H_wave_SF = fft(avg_wave)./H_mic;
    result.l_unmod_SF_f2 = H_wave_SF(f2/10+1,1);
    result.l_unmodCM_f2 = H_wave_SF(f2/10+1,2);

    % get unmodulated DPOAE and SF0AE at DP primaries
    avg_wave = (wave(30800 + latency+1 : fs/10 + 30800 + latency,:) ...
        - wave(38000 + latency+1 : fs/10 + 38000 + latency,:))/2;
    H_wave_unmod_DPodd = fft(avg_wave)./H_mic;

    avg_wave = (wave(30800 + latency+1 : fs/10 + 30800 + latency,:) ...
        + wave(38000 + latency+1 : fs/10 + 38000 + latency,:))/2;
    H_wave_unmod_DPeven = fft(avg_wave)./H_mic;

    result.l_unmod_2f1_f2 = H_wave_unmod_DPodd((2*f1-f2)/10+1,1);
    result.l_unmodCM_2f1_f2 = H_wave_unmod_DPodd((2*f1-f2)/10+1,2);
    result.l_unmod_f2_f1 = H_wave_unmod_DPodd((f2-f1)/10+1,1);
    result.l_unmodCM_f2_f1 = H_wave_unmod_DPodd((f2-f1)/10+1,2);
    result.l_unmod_2f2_f1 = H_wave_unmod_DPodd((2*f2-f1)/10+1,1);
    result.l_unmodCM_2f2_f1 = H_wave_unmod_DPodd((2*f2-f1)/10+1,2);
    result.l_unmod_SF_DPf1 = H_wave_unmod_DPodd(f1/10+1,1)-result.suppr_f1;
    result.l_unmodCM_DPf1 = H_wave_unmod_DPodd(f1/10+1,2);
    result.l_unmod_SF_DPf2 = H_wave_unmod_DPeven(f2/10+1,1)-result.suppr_f2;
    result.l_unmodCM_DPf2 = H_wave_unmod_DPeven(f2/10+1,2);
    result.l_unmod_2f1 = H_wave_unmod_DPeven(2*f1/10+1,1);
    result.l_unmodCM_2f1 = H_wave_unmod_DPeven(2*f1/10+1,2);

    %% biasing tone ON %%
    %% simulataneously presented primaries (modDP)
    modDP1_start = 74000 + latency;
    avg_wave = zeros(fs/10,2);
    for q = 0:7
        avg_wave = avg_wave + wave(q*6*fs/10+modDP1_start+1 : (q*6+1)*fs/10+modDP1_start,:);
    end
    avg_wave = avg_wave/8;
    result.H_modDP_1=fft(avg_wave)./H_mic;

    modDP2_start = 81200 + latency;
    avg_wave = zeros(fs/10,2);
    for q = 0:7
        avg_wave = avg_wave + wave(q*6*fs/10+modDP2_start+1 : (q*6+1)*fs/10+modDP2_start,:);
    end
    avg_wave = avg_wave/8;
    result.H_modDP_2=fft(avg_wave)./H_mic;

    H_mod_DPodd = (result.H_modDP_1-result.H_modDP_2)/2;
    H_mod_DPeven = (result.H_modDP_1+result.H_modDP_2)/2;

    % BT time course
    result.l_BT_ec = result.H_modDP_2(f_BT/10+1);
    H_f_BT = zeros(l,1);
    H_f_BT(f_BT/10+1) = result.l_BT_ec;
    course = real(ifft(H_f_BT));
    result.course_BT = course(1:fs/f_BT)/max(course);

    % 2f1-f2 suppression time courses
    mod_mask = zeros(l,1);
    mod_mask(mod_2f1_f2) = 1;

    result.H_mod_2f1_f2 = H_mod_DPeven(mod_2f1_f2,1);
    course = l * ifft(H_mod_DPeven(:,1) .* mod_mask);
    result.mod_2f1_f2_course = course(1:fs/f_BT);

    result.H_modCM_2f1_f2 = H_mod_DPeven(mod_2f1_f2,2);
    course = l * ifft(H_mod_DPeven(:,2) .* mod_mask);
    result.modCM_2f1_f2_course = course(1:fs/f_BT);

    % f2-f1 suppression time courses
    mod_mask = zeros(l,1);
    mod_mask(mod_f2_f1) = 1;

    result.H_mod_f2_f1 = H_mod_DPodd(mod_f2_f1,1);
    course = l * ifft(H_mod_DPodd(:,1) .* mod_mask);
    result.mod_f2_f1_course = course(1:fs/f_BT);

    result.H_modCM_f2_f1 = H_mod_DPodd(mod_f2_f1,2);
    course = l * ifft(H_mod_DPodd(:,2) .* mod_mask);
    result.modCM_f2_f1_course = course(1:fs/f_BT);

    % 2f2-f1 suppression time courses
    mod_mask = zeros(l,1);
    mod_mask(mod_2f2_f1) = 1;

    result.H_modCM_2f2_f1 = H_mod_DPodd(mod_2f2_f1,2);
    course = l * ifft(H_mod_DPodd(:,2) .* mod_mask);
    result.modCM_2f2_f1_course = course(1:fs/f_BT);

    result.H_mod_2f2_f1 = H_mod_DPodd(mod_2f2_f1,1);
    course = l * ifft(H_mod_DPodd(:,1) .* mod_mask);
    result.mod_2f2_f1_course = course(1:fs/f_BT);

    % 2f1 suppression time courses
    mod_mask = zeros(l,1);
    mod_mask(mod_2f1) = 1;

    result.H_modCM_2f1 = H_mod_DPeven(mod_2f1,2);
    course = l * ifft(H_mod_DPeven(:,2) .* mod_mask);
    result.modCM_2f1_course = course(1:fs/f_BT);

    result.H_mod_2f1 = H_mod_DPeven(mod_2f1,1);
    course = l * ifft(H_mod_DPeven(:,1) .* mod_mask);
    result.mod_2f1_course = course(1:fs/f_BT);

    % time courses of f1-SFOAE produced during simulataneously presented primaries and biasing tone
    result.l_SF_modDPf1 = H_mod_DPodd(f1/10+1,1) - result.suppr_f1;
    mod_mask = zeros(l,1);
    mod_mask(mod_f1) = 1;

    result.H_modCM_DPf1 = H_mod_DPodd(mod_f1,2);
    course = l * ifft(H_mod_DPodd(:,2) .* mod_mask);
    result.modCM_DPf1_course = course(1:fs/f_BT);

    H = H_mod_DPodd(:,1) .* mod_mask;
    H(f1/10+1)= result.l_SF_modDPf1;
    result.H_mod_DPf1 = H(mod_f1);
    course = l * ifft(H);
    result.modSF_DPf1_course = course(1:fs/f_BT);

    % time courses of f2-SFOAE produced during simulataneously presented primaries and biasing tone
    result.l_SF_modDPf2 = H_mod_DPeven(f2/10+1,1) - result.suppr_f2;
    mod_mask = zeros(l,1);
    mod_mask(mod_f2) = 1;

    result.H_modCM_DPf2 = H_mod_DPeven(mod_f2,2);
    course = l * ifft(H_mod_DPeven(:,2) .* mod_mask);
    result.modCM_DPf2_course = course(1:fs/f_BT);

    H = H_mod_DPeven(:,1) .* mod_mask;
    H(f2/10+1)= result.l_SF_modDPf2;
    result.H_mod_DPf2 = H(mod_f2);
    course = l * ifft(H);
    result.modSF_DPf2_course = course(1:fs/f_BT);


    %% only f2 primary presented (modSF)
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
    result.H_mod_SF = (H_modSF_1+H_modSF_2)/2;
    result.l_modSF_f2 = result.H_mod_SF(f2/10+1) - result.suppr_f2;

    % SFOAE time course
    mod_mask = zeros(l,1);
    mod_mask(mod_f2) = 1;

    result.H_modCM_f2 = result.H_mod_SF(mod_f2,2);
    course = l * ifft(result.H_mod_SF(:,2) .* mod_mask);
    result.modCM_f2_course = course(1:fs/f_BT);

    H = result.H_mod_SF(:,1) .* mod_mask;
    H(f2/10+1)= result.l_modSF_f2;
    result.H_mod_f2 = H(mod_f2);
    course = l * ifft(H);
    result.modSF_f2_course = course(1:fs/f_BT);
end


%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~plot_flag, return, end
plot_prmSet(results, prmSet,result,1,scrsz(monitor,:))

%plot time series of relevant spectral line levels
h_series = figure(get_figure_h([scrsz(monitor,1)+round(scrsz(monitor,3)/4*3),scrsz(monitor,2)+90,...
    round(scrsz(monitor,3)/4),round(scrsz(monitor,4)-150)])); clf
set(h_series,'Name',['Series' results.header.title])

if first==1
    clf
    series.H_mod_2f1_f2  = [];
    series.H_mod_f2_f1 = [];
    series.H_modCM_f2_f1 =[];
    series.H_modCM_2f1_f2 = [];
    series.H_mod_f2 = [];
    series.l_unmod_SF_f2 = [];
    series.H_modCM_f2 =[];
    series.l_unmodCM_f2 =[];
    series.CAP_ampl =[];
else
    series = getappdata(h_series,'series');
end

series.H_mod_2f1_f2 = [series.H_mod_2f1_f2 result.H_mod_2f1_f2];
series.H_modCM_2f1_f2 = [series.H_modCM_2f1_f2 result.H_modCM_2f1_f2];
series.H_mod_f2_f1 = [series.H_mod_f2_f1 result.H_mod_f2_f1];
series.H_modCM_f2_f1 = [series.H_modCM_f2_f1 result.H_modCM_f2_f1];
series.H_mod_f2 = [series.H_mod_f2 result.H_mod_f2];
series.l_unmod_SF_f2 = [series.l_unmod_SF_f2 result.l_unmod_SF_f2];
series.H_modCM_f2 = [series.H_modCM_f2 result.H_modCM_f2];
series.l_unmodCM_f2 = [series.l_unmodCM_f2 result.l_unmodCM_f2];
series.CAP_ampl = [series.CAP_ampl result.CAP_ampl];

setappdata(h_series,'series',series)
plot_time_series_spectral_lines(series)








