%% RCA analysis demo script - frequency domain sweep data

% This script explains how to run RCA, first on an individual subject level
% and then on the group level.
% For it to work, you will need to download / clone the SVNDL git
% repositories 'rcaBase' and 'mrC', and you can find them here:

% https://github.com/svndl/rcaBase
% https://github.com/svndl/mrC

% You will also need to have exported your EEG data in PowerDiva. This
% script assumes these data will be in .txt files and in RLS format (it can
% also handle DFT).

% This script is based around the functions 'rcaSweep' and 'aggregateData'
% to analyse SSVEP sweep data in the frequency domain.

% I have structured my analysis directory with the following subdirs:
% > code        contains rca toolboxes (don't need them here - can also be 
%               in git dir) and analysis scripts
% > exports     contains subdirs for each participant, each containing the 
%               PowerDiva exports
% > figures     for figures generated in Plot_RCA_Data.m
% > outputs     for outputs automatically generated in the RCA analysis

%% Clean up and add toolboxes to Matlab search path 

% Housekeeping
clear all; close all; clc;

% Add the RCA toolboxes
if ~exist('code_folder','var')
    code_folder = '/Users/kohler/code/';
    addpath(genpath(sprintf('%s/git/export_fig',code_folder)),'-end');
    addpath(genpath(sprintf('%s/git/rcaBase',code_folder)),'-end');
    addpath(genpath(sprintf('%s/git/mrC',code_folder)),'-end');
    addpath(genpath(sprintf('%s/git/schlegel/matlab_lib/misc',code_folder)),'-end');
    addpath(genpath(sprintf('%s/git/schlegel/matlab_lib/figure',code_folder)),'-end');
else
end


%% Set up some analysis parameters and variables

print_figures  = true;                                                                             % set to true if you want to automatically print the rca topo figures and save them into the'outputs' directory
save_name  = 'rcaData_Lochy_Cond3';    % name of output                                            % set to the desired ouput file name - this will contain one big data struct with the results from rcaSweep

data_location = '/Volumes/Seagate Backup Plus Drive/2019_synapse_data'; % change to your data directory
out_location = '/Volumes/Seagate Backup Plus Drive/2019_synapse_data/outputs';     % set to your output directory
fig_location = '/Volumes/Seagate Backup Plus Drive/2019_synapse_data/figures';   

% This loop will find all the the exported data folders in your data
% directory. Note I assume that all relevant folders start with the nl-
% identifier; you may need to change this depending on your own naming
% conventions.
% rcaSweep takes in a cell array specifying the exact locations of each
% participants' dataset. So, we need to write out a cell array that
% contains these.

folder_names = subfolders(sprintf('%s/BLC*',data_location),1); % get full path of all sub folders
sub_count = 0; 
for f = 1:length(folder_names)
    cur_folder = subfolders(sprintf('%s/*_EEGSsn', folder_names{f}),1);
    if cur_folder{1} ~= 0
        sub_count = sub_count + 1;
        new_names{sub_count} = sprintf('%s/Exp_TEXT_HCN_128_Avg', cur_folder{1});
    else
    end
end
folder_names = new_names; clear new_names;

% variable 'pathnames' should now contain all of the relevant directories

% Now we can set some more analysis-specific variables. For this part it is
% important to know how you want to analyse your data - i.e. how do you
% want the RCA filters to be constructed? Do you want to include all
% conditions, and all harmonics to analyse?

% this analysis is structured for Sweep analysis

binsToUse     = 1:10;          % Sweep indices of bins to include in analysis (the values must be present in the bin column of all DFT/RLS exports)
% binsToUse     = 0;           % Would typically do this though indices of bins to include in analysis (the values must be present in the bin column of all DFT/RLS exports)
freqsToUse    = [1 2];         % indices of frequencies to include in analysis (the values must be present in the frequency column of all DFT/RLS exports)... these freq were set when importing eegssn
                               % Chose freq 1 and 2, because 3Hz is carrier, so you want nf1 clean
condsToUse    = [3];           % if you want to include all conditions, create a vector here listing all condition numbers
                               % only look at condition 1 first, if you want to look at all 3 conditions
                               % together, do [1 2 3]
nReg          = 7;                          % RCA regularization constant (7-9 are typical values, but see within-trial eigenvalue plot in rca output); should always be bigger than the number of components. 
%trialsToUse   = [:];                       % which trials do you want to use? I've commented here so that it defaults to all trials it can find for each participant.
nComp         = 6;                          % number of RCs that you want to look at (3-5 are good values, but see across-trial eigenvalue plot in rca output)
chanToCompare = 65;                         % channel to use for a performance evaluation, we typically use 75 but should be determined on the basis of a representative 'good electrode' for your study - we use this just as a sanity check
dataType      = 'RLS';                      % can also be 'DFT' if you have DFT exports

rcPlotStyle = 'matchMaxSignsToRc1';         % not req'd. see 'help rcaRun', can be: 'matchMaxSignsToRc1' (default) or 'orig'

forceSourceData = 0;                        % set to true if you have source space; we didn't do this in our study

%% Run individual level RCA
% Hopefully, you won't need to change anything in this section

save_path = [out_location,'/',save_name]; % returns the full file path to your output folder

% Run some checks - have you analysed before / does the output file already
% exist?
launchAnalysis = false;
if isempty(dir([save_path,'.mat'])) % if .mat file specified by saveFileNamePath does not exist, go ahead and run analysis
   launchAnalysis = true; 
end

%% Run on the first 3 conditions (all oddball)
% (1) pseudo-fonts with word oddball
% (2) pseudo-fonts with non-word oddball 
% (3) non-words with word oddball

use_bins = 0;          % indices of bins to include in analysis (the values must be present in the bin column of all DFT/RLS exports)
%use_bins = 1:10;      % would typically do this though indices of bins to include in analysis (the values must be present in the bin column of all DFT/RLS exports)
use_freqs = 1:2;       % indices of harmonics to include in analysis 
                       % these will index the harmonics present in the DFT/RLS export
                       % Because 3Hz is the carrier, and 1 Hz is the
                       % oddball, we are including only the first two
                       % harmonics of the oddball
                        
use_conds = 1;         % if you want to include all conditions, create a vector here listing all condition numbers
                       % only look at condition 1 first, if you want to look at all 3 conditions
use_trials = [];

use_subs = [];
                       
n_reg = 7;             % RCA regularization constant (7-9 are typical values, but see within-trial eigenvalue plot in rca output); should always be bigger than the number of components. 
n_comp = 6;            % number of RCs that you want to look at (3-5 are good values, but see across-trial eigenvalue plot in rca output)
comp_channel = 65;     % channel to use for a performance evaluation, we typically use 75 but should be determined on the basis of a representative 'good electrode' for your study - we use this just as a sanity check
data_type = 'RLS';     % can also be 'DFT' if you have DFT exports

force_reload = false;  % if true, reload the data text files and generate
                       % new .mat files for quick loading of the data
                       % set to true if you have
                       % re-exported the data from xDiva    
rca_odd = rcaSweep(folder_names, use_bins, use_freqs, use_conds, use_trials, use_subs, n_reg, n_comp, data_type, comp_channel, force_reload);

if print_figures, print('-dpsc','~/Desktop/rca_odd.ps','-append'), end

% THIS SECTION OF THE CODE: 
% returns the structure 'rca_odd', whose dimensions are like this: 
% There will be as many 'structs' as HARMONICS you have analysed in the
% input data
% The output in W and A will be nChannels (128) x n RCA components

%% DO AVERAGE AND STATISTICS ON RCA STRUCT
% Now that we have run our analyses on individual subjects, we can generate
% group-level statistics by looping across these outputs using aggregateData.

keep_conditions = true;      % if false, will average across conditions
errorType      = 'SEM';      % can also be 'none' or '95CI' for 95% confidence intervals
trialError     = false;      % if true, will compute trial-wise errors and statistics
doNR           = [];         % Naka-Rushton fitting - turned off for now, but takes a logical array structured like harmonic x rc x condition

% loop over harmonics and aggregate data!
% note that the input and output is the same structure
% you are simply modifying the structure, adding additional fields
for f = 1:length(rca_odd)    
    fprintf('\n ... Harmonic no. %0.0f ...\n',f);
    rca_odd(f) = aggregateData(rca_odd(f), keep_conditions, errorType, trialError, doNR);
    % called like this:
    % aggregateData(rca_struct, keep_conditions, error_type, trial_wise, do_nr)
    % more info in 'help'
end

%% INSPECT THE RCA STRUCT 
h_idx = 1;
if strcmp(rca_odd(h_idx).settings.rcaSubs, 'all')
    num_subs = 'all';
else
    num_subs = num2str(length(rca_odd(h_idx).settings.rcaSubs));
end

if strcmp(rca_odd(h_idx).settings.rcaConds, 'all')
    num_conds = sprintf('all %d', num2str(size(rca_odd(h_idx).data, 1)));
else
    num_conds = num2str(length(rca_odd(h_idx).settings.rcaConds));
end

fprintf('\n I am an rca struct, created on %s, using %s data from %s subjects, and %s conditions. \n', ...
    rca_odd(h_idx).settings.runDate, ...
    rca_odd(h_idx).settings.data_type, ...
    num_subs, ...
    num_conds);

rca_freqs = rca_odd(h_idx).settings.freqLabels(rca_odd(h_idx).settings.rcaFreqs);
cur_freq = rca_odd(h_idx).settings.freqLabels{rca_odd(h_idx).settings.freqIndices};
str_to_print = ' I was made using data from';
for r = 1:(length(rca_freqs)-1)
    str_to_print = [str_to_print, sprintf(' %s,', rca_freqs{r}) ];
end
str_to_print = [str_to_print, sprintf(' and %s. This part of the array of structs contain data from %s.\n', rca_freqs{end}, cur_freq) ];
fprintf('%s', str_to_print)

%% FIGURE PARAMS
l_width = 1;
f_size = 12;
text_opts = {'fontweight','normal','fontname','Helvetica','fontsize', f_size};
gca_opts = [{'tickdir','out','ticklength',[0.0200,0.0200],'box','off','linewidth',l_width}, text_opts];

%% PLOT ERROR ELLIPSE
close all;
figure
% Specify which bin to plot (usually the average)
bin_to_plot = 'ave';
% Specify which RC component to plot
rc_idx = 1;
% Specify which harmonic to plot (always 1)
h_idx = 1;

c_idx = 1;

% Concatenate the real signal and the imaginary signal together 
% (bin x harmonic x subject x condition)
xy_vals = cat(2, squeeze(rca_odd(h_idx).subjects.real_signal(end, : , rc_idx, :, c_idx)), ...
                 squeeze(rca_odd(h_idx).subjects.imag_signal(end, :,  rc_idx, :, c_idx)));
% Set na vals so they'll be ignored
nan_vals = sum(isnan(xy_vals),2) > 0;

% Compute ellipse based on standard deviation and standard error of the
% mean - SEM is usually smaller and STD is usually bigger
[~,~,~,std_ellipse] = fitErrorEllipse(xy_vals(~nan_vals,:),'1STD');
[~,~,~,sem_ellipse] = fitErrorEllipse(xy_vals(~nan_vals,:),'SEM');

% Project xy_vals into vector space
vector_means = nanmean(xy_vals);
hold on
% Plot ellipses
plot(sem_ellipse(:,1), sem_ellipse(:,2),'-','linewidth',l_width,'Color', 'r');
plot(std_ellipse(:,1), std_ellipse(:,2),'-','linewidth',l_width,'Color', 'r');

% Plot individual subject data 
plot(xy_vals(~nan_vals,1),xy_vals(~nan_vals,2),'sq','linewidth',l_width,'Color', 'g');

% Plot vector amplitude
p_h = plot([0 vector_means(1)],[0 vector_means(2)],'-','linewidth',l_width,'Color','b');

% Create grid based on scatter of data
x_max = ceil(max(abs(xy_vals(~nan_vals,1))));
y_max = ceil(max(abs(xy_vals(~nan_vals,2))));
xlim([-x_max, x_max])
ylim([-y_max, y_max])
axis square 
plot([-x_max, x_max], zeros(1,2), 'k-', 'linewidth',l_width)
plot(zeros(1,2), [-y_max, y_max], 'k-', 'linewidth',l_width)
set(gca, gca_opts{:});
%% %% %% FLIP and REORDER RCA TOPOGRAPHIES and PHASE DATA
close all;
rc_idx = 1;
h_idx = 1;
c_idx = 1;
% Flip the RCA components and reorder them so component 1 becomes component
% 2. rcaOut, newOrder, index of components to flip
rca_new = rca_odd;

for h = 1:length(use_freqs)
    rca_new(h) = flipSwapRCA(rca_odd(h), [], rc_idx);
    rca_new(h) = aggregateData(rca_new(h), keep_conditions, errorType, trialError, doNR);
end
%% PLOT IT!
% Ellipse parameters
ellipseType = ['1STD' 'SEM'];

figure;
set(gcf, 'units', 'centimeters'); 
fig_x = 40; fig_y = 40;
fig_pos = get(gcf, 'position');
fig_pos(3:4) = [fig_x, fig_y];
fig_pos = set(gcf, 'position', fig_pos);

num_rows = 3;
num_cols = 3;
x_max = 30;
y_max = 10;
for z = 1:num_rows
    if z == 1  
        cur_rca = rca_odd;
        cur_idx = rc_idx;
        cur_w = rca_odd(h_idx).A(:,cur_idx);
    elseif z == 2
        cur_rca = rca_new;
        cur_idx = rc_idx;
        cur_w = rca_new(h_idx).A(:,cur_idx);
    else
        cur_rca = rca_new;
        cur_idx = n_comp + 1; % comparison channel
        n_electrodes = size(rca_odd(h_idx).W, 1);
        cur_w = zeros(n_electrodes,1); cur_w(chanToCompare) = 1; 
        x_max = 5;
        y_max = 2;
    end
    if z < 3
        % plot topographies
        subplot(num_rows, num_cols, 1+(z-1)*num_cols); hold on;
        mrC.plotOnEgi(cur_rca(h_idx).A(:,cur_idx))
        set(gca, gca_opts{:});
        hold off
    else
        subplot(num_rows, num_cols, 1+(z-1)*num_cols); hold on;
        mrC.plotOnEgi(ones(n_electrodes, 1)*.5, [], [], chanToCompare, true);
        colormap(gca, 'gray');
        set(gca, gca_opts{:});
        hold off
    end
    % plot frequency domain ellipse plots
    subplot(num_rows, num_cols, 2+(z-1)*num_cols)
    hold on
    xy_vals = cat(2, squeeze(cur_rca(h_idx).subjects.real_signal(end, : , cur_idx, :, c_idx)), ...
        squeeze(cur_rca(h_idx).subjects.imag_signal(end, :,  cur_idx, :, c_idx)));
    nan_vals = sum(isnan(xy_vals),2) > 0;
    % compute ellipse based on standard deviation and standard error of the
    % mean - SEM is usually smaller and STD is usually bigger
    [~,~,~,std_ellipse] = fitErrorEllipse(xy_vals(~nan_vals,:),'1STD');
    [~,~,~,sem_ellipse] = fitErrorEllipse(xy_vals(~nan_vals,:),'SEM');

    % Project xy_vals into vector space
    vector_means = nanmean(xy_vals);
    % Plot ellipses
    plot(sem_ellipse(:,1), sem_ellipse(:,2),'-','linewidth',l_width,'Color', 'r');  
    plot(std_ellipse(:,1), std_ellipse(:,2),'-','linewidth',l_width,'Color', 'r');

    % Plot individual subject data 
    plot(xy_vals(~nan_vals,1),xy_vals(~nan_vals,2),'sq','linewidth',l_width,'Color', 'g');

    % Plot vector amplitude
    p_h = plot([0 vector_means(1)],[0 vector_means(2)],'-','linewidth',l_width,'Color','b');

    % Create grid based on scatter of data
    % x_max = ceil(max(abs(xy_vals(~nan_vals,1))));
    % y_max = ceil(max(abs(xy_vals(~nan_vals,2))));
    
    xlim([-x_max, x_max])
    ylim([-y_max, y_max])
    axis square 
    plot([-x_max, x_max], zeros(1,2), 'k-', 'linewidth',l_width)
    plot(zeros(1,2), [-y_max, y_max], 'k-', 'linewidth',l_width)
    set(gca, gca_opts{:});
    hold off

    %Plot time domain
    subplot(num_cols, num_rows, 3+(z-1)*num_cols)
    hold on
    axx_rca_full = axxRCAmake(folder_names, cur_w);
    x_vals = linspace(0,1,420+1);
    x_vals = x_vals(2:end); 

    % Define colors (red   green   blue)
    colors = [1,0,0; 0,1,0; 0,0,1; 0 1 1];

    axx_mean = nanmean(cell2mat(axx_rca_full.Projected(c_idx,:)),2);
    axx_stderr = nanstd(cell2mat(axx_rca_full.Projected(c_idx,:)),0,2)./sqrt(length(folder_names));
    ph(c_idx) = plot(x_vals(1:length(axx_mean)), axx_mean, '-', 'linewidth', 2, 'color', colors(c_idx,:));
    ty_max = 25;
    ylim([-ty_max, ty_max])
    hold on; 
    ErrorBars(x_vals(1:length(axx_mean))', axx_mean, axx_stderr,'color',colors(c_idx,:));
    set(gca, gca_opts{:});
end

    
%%



    


subplot(2, 2, 4)
% Concatenate the real signal and the imaginary signal together 
% (bin x harmonic x subject x condition)
xy_new = cat(2, squeeze(rca_new(h).subjects.real_signal(end, : , rc_idx, :)), ...
                 squeeze(rca_new(h).subjects.imag_signal(end, :,  rc_idx, :)));

nan_new = sum(isnan(xy_new),2) > 0;

% Compute ellipse based on standard deviation and standard error of the
% mean - SEM is usually smaller and STD is usually bigger
[~,~,~,std_ellipse] = fitErrorEllipse(xy_new(~nan_new,:),'1STD');
[~,~,~,sem_ellipse] = fitErrorEllipse(xy_new(~nan_new,:),'SEM'); 


% Project xy_vals into vector space
new_vector_means = nanmean(xy_new);
hold on
% Plot ellipses
plot(sem_ellipse(:,1), sem_ellipse(:,2),'-','linewidth',l_width,'Color', 'r');
plot(std_ellipse(:,1), std_ellipse(:,2),'-','linewidth',l_width,'Color', 'r');

% Plot individual subject data 
plot(xy_new(~nan_new,1),xy_new(~nan_new,2),'sq','linewidth',l_width,'Color', 'g');

% Plot vector amplitude
p_h_new = plot([0 new_vector_means(1)],[0 new_vector_means(2)],'-','linewidth',l_width,'Color','b');

% Create grid based on scatter of data
new_x_max = ceil(max(abs(xy_new(~nan_new,1))));
new_y_max = ceil(max(abs(xy_new(~nan_new,2))));
xlim([-15, 15])
ylim([-8, 8])
axis square 
plot([-15, 15], zeros(1,2), 'k-', 'linewidth',l_width)
plot(zeros(1,2), [-8, 8], 'k-', 'linewidth',l_width)


%% Pull out phase significance - another way 

% Shows average t_sig for each component for the first harmonic
squeeze(rcaStruct(1).stats.t_sig(end, 1, :));

% Shows average t2_sig for each component for the first harmonic
squeeze(rcaStruct(1).stats.t2_sig(end, 1, :));

% Can also look at this bin-by-bin
squeeze(rcaStruct(1).stats.t2_sig(:, 1, :))

%% LOAD IN AXX DATA AND PLOT THEM FOR AN FREQ DOMAIN-DEFINED RC
close all;
rc_idx = 2;
h_idx = 1;
cur_w = rca_new(h_idx).A(:,rc_idx);
axx_rca_full = axxRCAmake(folder_names, cur_w);

x_vals = linspace(0,1,421);
x_vals = x_vals(2:end); 

% red   green   blue
colors = [1,0,0; 0,1,0; 0,0,1];

figure; clf
for c = 1:(size(axx_rca_full.Projected,1)-1) % All conditions except Vernier
    axx_mean = nanmean(cell2mat(axx_rca_full.Projected(c,:)),2);
    axx_stderr = nanstd(cell2mat(axx_rca_full.Projected(c,:)),0,2)./sqrt(length(folder_names));
    ph(c) = plot(x_vals, axx_mean, '-', 'linewidth', 2, 'color', colors(c,:));
    hold on; 
    ErrorBars(x_vals', axx_mean, axx_stderr,'color',colors(c,:));
end

legend(ph, 'Condition 1', 'Condition 2', 'Condition 3');

%% Component amplitude plot - in progress

comp_n = [1 2 3 4 5 6]


for h  = 1:length(use_freqs)
    for c = 1:length(nComp)
    % Find the mean amplitude of each component 
    curr_amp = rcaStruct(h).mean.amp_signal(end, : , c); %currently just for 
    curr_errlo = rcaStruct(h).stats.amp_lo_err(end, : , c);
    curr_errhi = rcaStruct(h).stats.amp_up_err(end, : , c);
    end
    plot(comp_n, curr_amp);
    
end
    
  

      
