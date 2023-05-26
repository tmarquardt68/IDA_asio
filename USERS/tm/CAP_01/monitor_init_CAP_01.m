function monitor_settings = monitor_init(results)

fs = results.stimulus.sample_rate;
monitor_settings.bp_filter = designfilt('bandpassiir','FilterOrder',20,...
    'PassbandFrequency1',350,'PassbandFrequency2',3000,'StopbandAttenuation1',60,...
    'PassbandRipple',1,'StopbandAttenuation2',80,'SampleRate',fs);