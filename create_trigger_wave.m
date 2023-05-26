function spike_wave = create_trigger_wave(fs),
% waveform of simulated spike
% Part of IDA_asio Toolbox
% Copyright (C) Torsten Marquardt
% Terms of the GNU General Public License apply
% (www.http://www.fsf.org/licensing/licenses/gpl.html).


spike_wave = -sin([1:1e-3*fs]*1000/fs*2*pi).* ...
    (1 - cos([1:1e-3*fs]*1000/fs*2*pi))/1.5;
spike_wave = spike_wave./max(max(spike_wave),-min(spike_wave));