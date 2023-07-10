function CAP
freqs = logspace(log10(2000),log10(22627.41699796954/2),31);
tone_pips = generate_tone_pips(freqs,10);
plot(tone_pips(1:end))

function tone_pips = generate_tone_pips(freqs,rep_rate)
global SAMPLE_RATE 
load()
cos_ramp_cycles = 3; % cycles
total_duration = 1/rep_rate*SAMPLE_RATE;
tone_pips = zeros(total_duration,numel(freqs));
for q = 1:numel(freqs)
    l_cos_ramp = round(cos_ramp_cycles/freqs(q) * SAMPLE_RATE);
    duration= 2*l_cos_ramp;
    tone_pips(:,q) = [...
        sin(linspace(0,freqs(q)/SAMPLE_RATE*duration*2*pi*(1-1/duration),duration))'... tone
         .* [cos(pi/2:(pi-pi/2)/(l_cos_ramp-1):pi).^2,cos(0:pi/2/(l_cos_ramp-1):pi/2).^2]';... enveleope
        zeros(total_duration-duration,1)]; % trailing silence
    Hdb_2 = 20*log10(abs(CURRENT_EAR(round(f2/10)+1,2)'));
end