function combi(varargin)
% combi_01('ARGNAME',ARGVAL,...) 
exp_description = 'Plays interleaved modSFOAE, modDPOAE, CM and Cap stimuli';
% options include:
% 'f2' : primary frequency (must be multiple of 20 Hz)
% 'l1' : dB above l2 (f1 is interanlly calucated to be f2/f1 ~ 1.2
% 'l2' : level of f2 
% 'l_suppr' : dB above l2 (to suppress SFOAE) (f_suppr is 40 Hz below f2)
% 'l_BT' : level of biasing tone
% 'n_rpts' : number of repeat (int)

%%
% PURPOSE:
% Stimuli generation for physiological experiments 
% v0.1 19-07-10 TM (Oldenburg)
% v1.1 18-10-12 TM (EI)
% v1.2 15-01-13 TM (EI)
% v2.0 21-05-13 GBC (EI)
% v2.2 13-10-15 TM (EI) 
% v3.0 8-03-22 TM (EI) monitor_combi_01
%=========================================================================
% First set baseline defaults.
n_rpts=1;
f2 = 7730;
f_BT = 60;

l1 = 10;
l2 = 65;
l_suppr = 10;
l_BT = [-40 75*ones(1,50) 90*ones(1,50) 75*ones(1,50 ) -40 ]+10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f2_old = f2;
f_BT_old = f_BT;
f2 = floor(f2/20)*20+10;
f_BT = floor(f_BT/20)*20;
f1 = floor(f2/1.25/20)*20;
f1_old = f1;
while mod(f2-f1,f_BT)==0
    f1 = f1+20;
end
if f2_old ~= f2
    warndlg(['f2 was ' num2str(f2_old) ' Hz;  Now: n*20+10  (' num2str(f2) ' Hz)'])
end
if f_BT_old ~= f_BT
    warndlg(['f_BT was ' num2str(f_BT_old) ' Hz; Now multiple of 20 (' num2str(f_BT) ' Hz)'])
end
if f1_old ~= f1
    warndlg(['F1 was not multiple of 20 Hz resulting f2-f1 being on n*20+10 grid. Now f1 = ' num2str(f1) ' Hz)'])
end
inter_onset_interval = []; % if empty: loop stimulus as fast as possible

% Now allow command line to overwrite any of the above.  In principle, this
% could allow you to accidentsally use the command line to change stuff 
% that's hidden in the config file. But if you are able to figure out how  
% to do that, you probably mean to do it!
warnopts(assignopts(who,varargin{:}));

cd(fileparts(which(mfilename))), % make sure it runs local stimulus file
load(['CONFIGS\' mfilename], 'stimuli_parameter');

%% BEGIN: DO NOT MODIFY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%% END: DO NOT MODIFY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stimulus.sample_rate = 48000;
stimulus.presentation_period = inter_onset_interval; % [s]

% parameter table definition
n=n+1;
parameters(n).stimulus_no = 1;
parameters(n).name = 'f_BT [Hz]';
parameters(n).unit = '';
parameters(n).values = f_BT; 
parameters(n).description = [];

n=n+1;
parameters(n).stimulus_no = 1;
parameters(n).name = 'f1 [Hz]';
parameters(n).unit = '';
parameters(n).values = f1; 
parameters(n).description = [];

n=n+1;
parameters(n).stimulus_no = 1;
parameters(n).name = 'f2 [Hz]';
parameters(n).unit = '';
parameters(n).values = f2; 
parameters(n).description = [];

n=n+1;
parameters(n).stimulus_no = 1;
parameters(n).name = 'l1 [+dB]';
parameters(n).unit = '';
parameters(n).values = l1;
parameters(n).description = [];

n=n+1;
parameters(n).stimulus_no = 1;
parameters(n).name = 'l2 [dB]';
parameters(n).unit = '';
parameters(n).values = l2;
parameters(n).description = [];

n=n+1;
parameters(n).stimulus_no = 1;
parameters(n).name = 'l_suppr [+dB]';
parameters(n).unit = '';
parameters(n).values = l_suppr;
parameters(n).description = [];

n=n+1;
parameters(n).stimulus_no = 1;
parameters(n).name = 'l_BT [dB]';
parameters(n).unit = '';
parameters(n).values = l_BT;
parameters(n).description = [];

% create parameter table
n=0;
for q2 = 1:size(parameters(2).values,2)
    for q3 = 1:size(parameters(3).values,2)
        for q4 = 1:size(parameters(4).values,2)
            for q5 = 1:size(parameters(5).values,2)
                for q6 = 1:size(parameters(6).values,2)
                    for q7 = 1:size(parameters(7).values,2)
                        for q8 = 1:size(parameters(8).values,2)
                            n = n+1;
                            parameter_table(n,:) = [n, ...
                                parameters(2).values(q2), ...
                                parameters(3).values(q3), ...
                                parameters(4).values(q4), ...
                                parameters(5).values(q5), ...
                                parameters(6).values(q6), ...
                                parameters(7).values(q7), ...
                                parameters(8).values(q8)];
                            original_table_order(n,:) = [n q2 q3 q4 q5];
                        end
                    end
                end
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
ida('run_experiment',mfilename,stimulus,parameter_table,exp_description)
%% END: DO NOT MODIFY

