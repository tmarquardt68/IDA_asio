function CAP(varargin)
% CAP_01('ARGNAME',ARGVAL,...) runs a frequency response area.  ARGNAME
% options include:
%
% 'freqs' : array of frequencies (Hz) to present (default: 200 Hz to 6400
%            Hz in 1/8 octave steps)
% 'levels' : array of levels (dB) to present (default: 40 to 90 in 10 dB steps)
% 'n_rpts' : number of repetitions (default 10)
% 'inter_onset_interval' : time (s) between stimulus presentations
%            (default: 0.4)
% 'animal_id' : ID of animal (string)
% 'track_no' : track number (int)
% 'pos_no' : position number (int)
% 'rpt_no' : repeat number (int)

%%
% PURPOSE:
% Stimuli generation for physiological experiments 
% v0.1 19-07-10 TM (Oldenburg)
% v1.1 18-10-12 TM (EI)
% v1.2 15-01-13 TM (EI)
% v2.0 21-05-13 GBC (EI)
% v2.2 13-10-15 TM (EI) 
% v3.0 8-03-22 TM (EI) 
%=========================================================================
exp_description = 'Checking required gaps between stimulis for full CAP recovery';

% First set baseline defaults.
n_rpts=1;
gaps = 50;
levels = 100;
tone_lengths = 0;
freqs = logspace(log10(2000),log10(22627.41699796954/2),31);

inter_onset_interval = 1;

%% BEGIN: DO NOT MODIFY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now allow command line to overwrite any of the above.  In principle, this
% could allow you to accidentsally use the command line to change stuff 
% that's hidden in the config file. But if you are able to figure out how  
% to do that, you probably mean to do it!
warnopts(assignopts(who,varargin{:}));

config_filename=mfilename;
cd(fileparts(which(config_filename))), % make sure it runs local stimulus file
load(config_filename, 'stimuli_parameter');

n=1;
parameters(n).stimulus_no = 0;
parameters(n).name = 'index';
parameters(n).unit = '';
parameters(n).description = 'Running index of stimulus table as generated below';
parameters(n).values = [];

parameter_table = [];
stimulus.raw = [];
stimulus.stimuli_parameters = stimuli_parameter;
stimulus.n_rpts = n_rpts;
stimulus.presentation_period = inter_onset_interval; % [s]
%% END: DO NOT MODIFY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stimulus.sample_rate = 48000;

% parameter table definition
n=n+1;
parameters(n).stimulus_no = 1;
parameters(n).name = 'freq [Hz]';
parameters(n).unit = '';
parameters(n).values = freqs; % need to round because freqs must be int
parameters(n).description = [];

n=n+1;
parameters(n).stimulus_no = 1;
parameters(n).name = 'tone_len [ms]';
parameters(n).unit = '';
parameters(n).values = tone_lengths; % need to round because freqs must be int
parameters(n).description = [];

n=n+1;
parameters(n).stimulus_no = 1;
parameters(n).name = 'level [dB]';
parameters(n).unit = '';
parameters(n).values = levels;
parameters(n).description = [];

n=n+1;
parameters(n).stimulus_no = 1;
parameters(n).name = 'gap [ms]';
parameters(n).unit = '';
parameters(n).values = gaps;
parameters(n).description = [];


% create parameter table
n=0;
for q2 = 1:size(parameters(2).values,2)
    for q3 = 1:size(parameters(3).values,2)
        for q4 = 1:size(parameters(4).values,2)
            for q5 = 1:size(parameters(5).values,2)
                n = n+1;
                parameter_table(n,:) = [n, ...
                    parameters(2).values(q2), ...
                    parameters(3).values(q3), ...
                    parameters(4).values(q4), ...
                    parameters(5).values(q5)];
                original_table_order(n,:) = [n q2 q3 q4 q5];
            end
        end
    end
end

% structure describing the scheme of generating presentation order (Still to be defined!)
stimulus.ordering_scheme.block.ordering_type = 0; % 1 - random

%% BEGIN: DO NOT MODIFY
stimulus.ordering_scheme.block.original_table_order = original_table_order; 
stimulus.parameters = parameters;
for q = 1:n_rpts
    if stimulus.ordering_scheme.block.ordering_type == 1
        interval_order(q,:) = randperm(size(parameter_table,1));
    else
        interval_order(q,:) = 1:size(parameter_table,1);
    end
end

stimulus.interval_order = interval_order;
% Start presentation
ida('run_experiment',config_filename,stimulus,parameter_table,exp_description)
%% END: DO NOT MODIFY

