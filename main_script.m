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
close all; clc;

% ! change the variables below to point to your git and data folder !
data_location = '/Volumes/BACKUP/Lochy_RCA';
git_folder = '~/code/git';

% Add the RCA toolboxes
if isempty(which('rcaRun'))
    addpath(genpath(fullfile(git_folder, 'export_fig')),'-end');
    addpath(genpath(fullfile(git_folder, 'rcaBase')),'-end');
    addpath(genpath(fullfile(git_folder, 'mrC')),'-end');
    if exist(genpath(fullfile(git_folder, 'matlab_lib')), 'dir')
        addpath(genpath(fullfile(git_folder, 'matlab_lib', 'misc')),'-end');
        addpath(genpath(fullfile(git_folder, 'matlab_lib', 'figure')),'-end');
    else
        % pjkohler's version
        addpath(genpath(fullfile(git_folder, 'schlegel', 'matlab_lib', 'misc')),'-end');
        addpath(genpath(fullfile(git_folder, 'schlegel', 'matlab_lib', 'figure')),'-end');
    end
else
end

%% Set up some analysis parameters and variables

print_figures  = true;                                                                             % set to true if you want to automatically print the rca topo figures and save them into the'outputs' directory
save_name  = 'rcaData_Lochy_Cond1';    % name of output                                            % set to the desired ouput file name - this will contain one big data struct with the results from rcaSweep

 % change to your data directory
out_location = fullfile(data_location, 'outputs');     % set to your output directory
fig_location = fullfile(data_location, 'figures');    

% This loop will find all the the exported data folders in your data
% directory. Note I assume that all relevant folders start with the nl-
% identifier; you may need to change this depending on your own naming
% conventions.
% rcaSweep takes in a cell array specifying the exact locations of each
% participants' dataset. So, we need to write out a cell array that
% contains these.

if ispc
    file_sep = '\';
else
    file_sep = '/';
end

folder_names = subfolders(sprintf('%s%snl*',data_location, file_sep), 1); % get full path of all sub folders

if ~folder_names{1}
    folder_names = subfolders(sprintf('%s%sblc*',data_location, file_sep), 1); % get full path of all sub folders
else
end

if ~folder_names{1}
    msg = 'No SSVEP data in data_location!';
    error(msg)
else
end

sub_count = 0;
for f = 1:length(folder_names)
    cur_folder = subfolders(sprintf('%s%s*_EEGSsn', folder_names{f}, file_sep), 1);
    if cur_folder{1} ~= 0
        sub_count = sub_count + 1;
        new_names{sub_count} = sprintf('%s%sExp_TEXT_HCN_128_Avg', cur_folder{1}, file_sep);
    else
    end
end
folder_names = new_names; clear new_names

% variable 'pathnames' should now contain all of the relevant directories

%% Run on the first 3 conditions (all oddball)
% (1) pseudo-fonts with word oddball
% (2) pseudo-fonts with non-word oddball 
% (3) non-words with word oddball

use_bins = 0;          % indices of bins to include in analysis (the values must be present in the bin column of all DFT/RLS exports)
%use_bins = [0, 2, 3, 7, 9];      % would typically do this though indices of bins to include in analysis (the values must be present in the bin column of all DFT/RLS exports)
use_freqs = [1:4, 6:9];       % indices of harmonics to include in analysis 
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
                       
return_all = true;    % if true (the default) returns data from of the data in the dataset,
                       % all bins, harmonics, conditions, subjects, trials,
                       % regardless of what subset of the data rca has been
                       % trained on. If false, only return data that rca
                       % was trained on. w
rca_test = rcaSweep...
    (folder_names, use_bins, use_freqs, use_conds, use_trials, use_subs, n_reg, n_comp, data_type, comp_channel, force_reload, return_all);

if print_figures, print('-dpsc', fullfile(data_location, 'figures', 'rca_test.ps'), '-append'), end

% THIS SECTION OF THE CODE: 
% returns the structure 'rca_test', whose dimensions are like this: 
% There will be as many 'structs' as HARMONICS you have analysed in the
% input data
% The output in W and A will be nChannels (128) x n RCA components

%% FIGURE PARAMETERS (will be used later)

l_width = 1;
f_size = 12;
text_opts = {'fontweight','normal','fontname','Helvetica','fontsize', f_size};
gca_opts = [{'tickdir','out','ticklength',[0.0200,0.0200],'box','off','linewidth',l_width, 'clipping', 'off'}, text_opts];

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
for f = 1:length(rca_test)    
    fprintf('\n ... Harmonic no. %0.0f ...\n',f);
    rca_test(f) = aggregateData(rca_test(f), keep_conditions, errorType, trialError, doNR);
    % called like this:
    % aggregateData(rca_struct, keep_conditions, error_type, trial_wise, do_nr)
    % more info in 'help'
end

%% INSPECT THE RCA STRUCT 
h_idx = 1;
if strcmp(rca_test(h_idx).settings.rcaSubs, 'all')
    num_subs = 'all';
else
    num_subs = num2str(length(rca_test(h_idx).settings.rcaSubs));
end

if strcmp(rca_test(h_idx).settings.rcaConds, 'all')
    num_conds = sprintf('all %d', num2str(size(rca_test(h_idx).data, 1)));
else
    num_conds = num2str(length(rca_test(h_idx).settings.rcaConds));
end

fprintf('\n I am an rca struct, created on %s, using %s data from %s subjects, and %s conditions. \n', ...
    rca_test(h_idx).settings.runDate, ...
    rca_test(h_idx).settings.data_type, ...
    num_subs, ...
    num_conds);

rca_freqs = rca_test(h_idx).settings.freqLabels(rca_test(h_idx).settings.rcaFreqs);
cur_freq = rca_test(h_idx).settings.freqLabels{rca_test(h_idx).settings.freqIndices};
str_to_print = ' I was made using data from';
for r = 1:(length(rca_freqs)-1)
    str_to_print = [str_to_print, sprintf(' %s,', rca_freqs{r}) ];
end
str_to_print = [str_to_print, sprintf(' and %s. This part of the array of structs contain data from %s.\n', rca_freqs{end}, cur_freq) ];
fprintf('%s', str_to_print)

%% LOAD IN AXX DATA
% set parameters
h_idx = 1;
cur_w = rca_test(h_idx).W; % note, h_idx does not matter here, because W is the same for all harmonics

% make axxRCA struct 
axx_test = rcaWaveProject(folder_names, cur_w);
% FANG: WRITE SOMETHING ABOUT HOW COOL rcaWaveProject is! 
% [dual functionality, input can be both pathnames or 
% a cell matrix of sensor data

%% ASSESS WHICH RC COMPONENTS SHOULD BE FLIPPED
% determine which components should be flipped based on the axx data
[flip_list, corr_list] = componentComparison(axx_test);
comp_to_flip = find(flip_list(1,:));
% FANG: WRITE SOME CODE TO EXPRESS THE CONSISTENTY OF COMPONENT SIGN ACROSS
% THE HIGH CORRELATIONS

%% FLIP and REORDER RCA TOPOGRAPHIES and PHASE DATA
close all;
% create rca_new struct to be flipped
rca_flipped = rca_test;

new_order = []; % note, you can use the function to flip your RC order as well
                % just provide the order here e.g. new_order = [2, 1, 3 ,4 5]

% flip rca components and data
for h = 1:length(rca_test)
    rca_flipped(h) = flipSwapRCA(rca_test(h), new_order, comp_to_flip);
    rca_flipped(h) = aggregateData(rca_flipped(h), keep_conditions, errorType, trialError, doNR);
end

cur_w = rca_flipped(h_idx).W; % note, h_idx does not matter here, because W is the same for all harmonics

% make axxRCA struct for flipped data
axx_flipped = rcaWaveProject(folder_names, cur_w);

%% DO TIME-DOMAIN STATISTICS USING AXX DATA
conds_to_compare = [1, 2];
% reorganize structure to take the difference of the waveforms 
diff_comp = cell2mat(permute(axx_flipped.rcaWave(conds_to_compare,:),[3,2,1]));
diff_comp = diff_comp(:,:,1) - diff_comp(:,:,2);

% set max number of permutations 
n_perms = 50000;

% choose type of test statistic - size, height, or mass - mass accounts for length & amplitude
test_stat = 'mass'; 
 
% run the ttest_permute_sstats function
[realH, realP, realT, corrH, critVal, supraThr, clustDistrib ]= ttest_permute_sstats(diff_comp, n_perms, test_stat);

%% ...  AND PLOT THEM FOR AN FREQ DOMAIN-DEFINED RC
close all;

% significance colors
p_colormap = jmaColors('pval');
p_colormap(end,:) = [1 1 1];

% set parameters
num_t = size(axx_flipped.rcaWave{conds_to_compare(1),1},1);
axx_xvals = linspace(0,1,420+1);
axx_xvals = axx_xvals(2:end);
axx_xvals = axx_xvals(1:num_t);
axx_xmin = 0; axx_xmax = max(axx_xvals);

% red   green   blue
colors = [1,0,0; 0,1,0; 0,0,1];

figure;
hold on;
% plot line at 0
plot([axx_xmin, axx_xmax], zeros(2,1), 'k-', 'linewidth', l_width)

axx_ymax = 0;

% plot averaged waveforms from each condition and standard errors
for c = 1:length(conds_to_compare) % All conditions except Vernier
    axx_mean = nanmean(cell2mat(axx_flipped.rcaWave(conds_to_compare(c),:)),2);
    axx_stderr = nanstd(cell2mat(axx_flipped.rcaWave(conds_to_compare(c),:)),0,2)./sqrt(length(folder_names));
    ph(c) = plot(axx_xvals, axx_mean, '-', 'linewidth', l_width, 'color', colors(c,:));
    ErrorBars(axx_xvals', axx_mean, axx_stderr,'color',colors(c,:));
    cur_max = ceil(max(abs([axx_mean-axx_stderr; axx_mean+axx_stderr])));
    if cur_max > axx_ymax
        axx_ymax = cur_max;
    else
    end
end

axx_ymin = - axx_ymax; axx_yunit = 1;
sig_pos(1) = axx_ymax;
sig_pos(2) = axx_ymax -( axx_ymax-axx_ymin) *.05; 

% plot corrected t-values
regionIdx = bwlabel(corrH);
for m=1:max(regionIdx)
    tmp = regionprops(regionIdx == m,'centroid');
    idx = round(tmp.Centroid(2));
    hTxt = text(axx_xvals(idx), sig_pos(1),'*','fontsize',36,'fontname','Helvetica','horizontalalignment','center','verticalalignment','cap');
end

% plot uncorrected t-values
cur_p = repmat( realP', 20,1 );
h_img = image([min(axx_xvals),max(axx_xvals)],[sig_pos(1),sig_pos(2)], cur_p, 'CDataMapping', 'scaled','Parent',gca);
colormap( gca, p_colormap );   
c_mapmax = .05+2*.05/(size(p_colormap,1));
set( gca, 'CLim', [ 0 c_mapmax ] ); % set range for color scale
set(gca, gca_opts{:});
uistack(h_img,'bottom')

% plot averaged difference of conditions
diff_mean = nanmean(diff_comp,2);
diff_h = plot(axx_xvals, diff_mean, '-', 'linewidth', l_width, 'color', colors(c+1,:));

% specify plot limits and legend
cond_labels = arrayfun(@(x) ['\itcondition ', num2str(x)], conds_to_compare, 'uni', false);
xlim([axx_xmin, axx_xmax]);
ylim([axx_ymin, axx_ymax]);
legend([ph, diff_h], [cond_labels, '\itdifference'], text_opts{:}, 'box', 'off', 'location', 'southwest');
set(gca, gca_opts{:}, 'ytick', axx_ymin:axx_yunit:axx_ymax, 'xtick', axx_xmin:.1:axx_xmax);
xlabel('{\ittime} (s)', text_opts{:});
ylabel('{\itamplitude} (\muV)', text_opts{:});

hold off
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
xy_vals = cat(2, squeeze(rca_test(h_idx).subjects.real_signal(end, : , rc_idx, :, c_idx)), ...
                 squeeze(rca_test(h_idx).subjects.imag_signal(end, :,  rc_idx, :, c_idx)));
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

%% COMPARE FLIPPED AND NON-FLIPPED DATA
% in the frequency and time domain
h_idx = 1;
c_idx = 1;
rc_idx = 1;
% ellipse parameters
ellipseType = ['1STD' 'SEM'];

% set figure parameters
figure;
set(gcf, 'units', 'centimeters'); 
fig_x = 40; fig_y = 40;
fig_pos = get(gcf, 'position');
fig_pos(3:4) = [fig_x, fig_y];
fig_pos = set(gcf, 'position', fig_pos);

num_rows = 3;
num_cols = 3;
x_max = 2;
y_max = 2;
ty_max = 20;
for z = 1:num_rows
    if z == 1  
        cur_rca = rca_test;
        cur_idx = rc_idx;
        cur_w = rca_test(h_idx).W(:,cur_idx);
    elseif z == 2
        cur_rca = rca_flipped;
        cur_idx = rc_idx;
        cur_w = rca_flipped(h_idx).W(:,cur_idx);
    else
        cur_rca = rca_test;
        cur_idx = n_comp + 1; % comparison channel
        n_electrodes = size(rca_test(h_idx).W, 1);
        cur_w = zeros(n_electrodes,1); cur_w(comp_channel) = 1; 
        x_max = 1;
        y_max = 1;
        ty_max = 4;
    end
    if z < 3
        % plot topographies
        subplot(num_rows, num_cols, 1+(z-1)*num_cols); hold on;
        mrC.plotOnEgi(cur_rca(h_idx).A(:,cur_idx));
        set(gca, gca_opts{:});
        hold off
    else
        subplot(num_rows, num_cols, 1+(z-1)*num_cols); hold on;
        mrC.plotOnEgi(ones(n_electrodes, 1)*.5, [], [], comp_channel, true);
%         cube = colorcube; 
%         cube = cube(250, :);
        colormap(gca, gray);
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
    axx_test_full = rcaWaveProject(folder_names, cur_w);
    x_vals = linspace(0,1,420+1);
    x_vals = x_vals(2:end); 

    % Define colors (red   green   blue)
    colors = [1,0,0; 0,1,0; 0,0,1; 0 1 1];

    axx_mean = nanmean(cell2mat(axx_test_full.rcaWave(c_idx,:)),2);
    axx_stderr = nanstd(cell2mat(axx_test_full.rcaWave(c_idx,:)),0,2)./sqrt(length(folder_names));
    ph(c_idx) = plot(x_vals(1:length(axx_mean)), axx_mean, '-', 'linewidth', 2, 'color', colors(c_idx,:));
    ylim([-ty_max, ty_max])
    hold on; 
    ErrorBars(x_vals(1:length(axx_mean))', axx_mean, axx_stderr,'color',colors(c_idx,:));
    set(gca, gca_opts{:});
end

%% Pull out phase significance - another way 

% Shows average t_sig for each component for the first harmonic
%squeeze(rca_test(1).stats.t_sig(end, 1, :));

% Shows average t2_sig for each component for the first harmonic
%squeeze(rcaStruct(1).stats.t2_sig(end, 1, :));

% Can also look at this bin-by-bin
%squeeze(rca_test(1).stats.t2_sig(:, 1, :))

%% COMPARE CONDITIONS IN THE FREQUENCY DOMAIN
% set parameters
rc_idx = 1;
bin_idx = 11;
h_idx = 1;
conds_to_compare = [1, 2];

% run student's t-test on projected amplitudes
project_data = squeeze(rca_flipped(h_idx).subjects.proj_amp_signal(bin_idx,:,rc_idx,:,conds_to_compare));
[project_T, project_P, ci, temp_stats] = ttest(project_data(:,1), project_data(:,2), 'dim', 1, 'tail', 'both');

% run Hotelling's T-squared test using real and imaginary values
real_data = squeeze(rca_flipped(h_idx).subjects.real_signal(bin_idx,:,rc_idx,:,conds_to_compare));
imag_data = squeeze(rca_flipped(h_idx).subjects.imag_signal(bin_idx,:,rc_idx,:,conds_to_compare));
combined_data = cat(3, real_data, imag_data);
combined_data = permute(combined_data, [1,3,2]);
t2results = tSquaredFourierCoefs(combined_data);

harm_text = rca_flipped(h_idx).settings.freqLabels{h_idx};
bin_text = rca_flipped(h_idx).settings.binLabels{bin_idx};

% create output of results
fprintf('\ncomparing conditions %d and %d, %s bin and %s harmonic\n', ...
    conds_to_compare(1), conds_to_compare(2), bin_text, harm_text);  
fprintf('-> student''s t-test: t(%d) = %.3f, p = %.3f\n', temp_stats.df, temp_stats.tstat, project_P);
fprintf('-> Hotelling''s t2: t2(%d, %d) = %.3f, p = %.3f \n\n', t2results.df1, t2results.df2, t2results.tSqrd, t2results.pVal)

%% Component amplitude plot - in progress

% comp_n = [1 2 3 4 5 6]
% 
% 
% for h  = 1:length(use_freqs)
%     for c = 1:length(nComp)
%     % Find the mean amplitude of each component 
%     curr_amp = rca_flipped(h).mean.amp_signal(end, : , c); %currently just for 
%     curr_errlo = rca_flipped(h).stats.amp_lo_err(end, : , c);
%     curr_errhi = rca_flipped(h).stats.amp_up_err(end, : , c);
%     end
%     bar(comp_n, curr_amp);
%     
% end
    
  

      
