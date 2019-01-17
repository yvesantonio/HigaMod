function acq = PA_import(acq_dir,show_plot,save_name)
%% PA_import - Import data from PoliArd acquisition
%
% Returns in a structure signal names, signals vectors, time vector
% from PoliArd folder acquisition.
% 
% acq = PA_import                                Wizard mode
% acq = PA_import(acq_dir)                       Acquisition folder
% acq = PA_import(acq_dir,show_plot)             Show plot (bool)
% acq = PA_import(acq_dir,show_plot,save_name)   Save to .MAT

% Oct 2015 - Damiano Milani - damiano.milani@polimi.it
% Politecnico di Milano - Department of Mechanical Engineering
% mecsys.mecc.polimi.org/mechlab


% Args parsing
choose_dir = false;
save_enable = true;

if nargin==0        % wizard mode
    choose_dir = true;
    show_plot = true;
    save_enable = true;
end

if nargin==1
    show_plot = false;
    save_enable = false;
end

if nargin==2
    save_enable = false;
end

%%
fprintf('\n IMPORT DATA FROM PoliArd ACQUISITION FILES \n')
fprintf('============================================\n\n')

if choose_dir
    % Choose directory from menu
    acq_dir = uigetdir('Select acquisition folder...');
    if ~acq_dir        %return if no file selected
        return
    end
end

file_dir = dir([acq_dir '\data*.txt']);
acq_list = {file_dir(:).name};


%% LOAD FILES
fprintf('* Acquisition #%s\n\n',acq_dir(end-14:end))

% Read info from config file
fprintf('* Opening config file...\n')
fileconf_dir = fullfile(acq_dir,'conf_rec.txt');
fileconf = fopen(fileconf_dir,'r');
C = textscan(fileconf, '%s', 'Delimiter', '\n');
fclose(fileconf);
head_target = '+ Signals names: ';
row_f = ~cellfun('isempty', strfind(C{1}, head_target));
names = strrep(C{1}{row_f},head_target,'');

acq.name = strsplit(names,', ');

fprintf([C{1}{1} '\n' C{1}{2} '\n# Signals: '])
fprintf('%s ',acq.name{:})

acq.data = [];
acq.time = [];
    
head_target = '+ Sampling frequency [Hz]: ';
row_f = ~cellfun('isempty', strfind(C{1}, head_target));
f_samp = str2double(strrep(C{1}{row_f},head_target,''));
fprintf('\n# Sampling frequency: %.2f\n\n',f_samp)

dt = 1/f_samp;

% Read info from acq files
for ff = 1:length(acq_list)
    fprintf('* Importing acq file %d/%d...\n',ff,length(acq_list))
    fileandpath = fullfile(acq_dir,acq_list{ff});
    
    %signal data
    datafile = importdata(fileandpath,'\t'); 
    acq.data = [acq.data; datafile];
    n=numel(acq.data(1,:));
    %time vector
    if ff==1
        t0 = 0;
    else
        t0 = time(end);
    end
    time = t0 + (dt*n:dt*n:dt*n*size(datafile,1));
    acq.time = [acq.time; time'];
end


%% PLOT DATA
if show_plot

    figure
    plot(acq.time,acq.data)
    grid on
    title('PoliArd Acquisition')
    legend(acq.name{:})
    xlabel('t [s]')
end


%% SAVE TO .MAT FILE
if save_enable
    if nargin <3  % *--*
        [save_file,save_path] = uiputfile('*.mat','Save data to MAT file',[acq_dir(end-14:end) '.mat']);
        if ~save_path        %return if no file selected
            return
        end
        save_name = fullfile(save_path,save_file);
    end
    
    save(save_name,'acq');
    fprintf('\n* Data saved in %s\n',save_name)
end
fprintf('============================================\n\n')