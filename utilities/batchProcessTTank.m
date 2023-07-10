function anCount=batchProcessTTank(vcFiles,varargin)

%% If no arguments passed, allow selection via graphical interface.
if nargin<1 || isempty(vcFiles)
    asTankPath=uigetdir('C:\TDT\OpenEx\MyProjects', ...
        'Select a tank:');
    vcBlockList=dir([asTankPath filesep 'Block*']);
    vcBlockList={vcBlockList.name};
    [~,asTankName]=fileparts(asTankPath);
    vcBlockChoices=listdlg('PromptString','Select blocks to process:',...
        'ListString',vcBlockList,'Name',asTankName,'OKString','Begin');
    if isempty(vcBlockChoices)
        return;
    end
    vcFiles=cell(1,length(vcBlockChoices));
    for n=1:length(vcBlockChoices)
        vcFiles{n}=[asTankPath filesep vcBlockList{vcBlockChoices(n)}];
    end
    asAnalogAnswer=questdlg('Process analogue data?',...
        'Analogue data','Process','Ignore','Process');
    switch asAnalogAnswer
        case 'Ignore'
            abDoAnalog=false;
        case 'Process'
            abDoAnalog=true;
    end
else
    if ~iscell(vcFiles)
        vcFiles={vcFiles};
    end
    abDoAnalog=true;
end

abDoPlotting=false;
warnopts(assignopts(who,varargin{:}));

[~,asTankName]=fileparts(asTankPath);
asPlotDir=['C:\Users\suser\Documents\MATLAB\IDA_asio_TDT\FIGS\' asTankName];
mkdir(asPlotDir)
anCount=0;
for f=1:length(vcFiles)
    fprintf(1,'Processing %s (%d of %d)...\n',vcFiles{f},f,length(vcFiles));
    [atResults,asOutFile]=processTTank(...
        vcFiles{f},'doSave',true,'doAnalog',abDoAnalog,varargin{:});
    if abDoPlotting
        vuHandles=eval(['plot_' atResults.header.config_filename '(atResults)']);
        [~,asOutFile]=fileparts(asOutFile);
        for n=1:length(vuHandles)
            print(vuHandles(n),'-dtiff','-r300',[asPlotDir filesep asOutFile '-' num2str(n)]);
            close(vuHandles(n));
        end
    end
end