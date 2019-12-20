%% add toolboxes to Matlab search path 
close all; clc;

% ! change the variables below to point to your git and data folder !
%data_location = '/Volumes/BACKUP/Lochy_RCA';
data_location = '/Volumes/GoogleDrive/My Drive/2019_synapse_data';
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

%% get indices for selecting sub-populations
reader_file = fullfile(data_location, 'BLC_MSbehavioral_12102019.csv');
reader_data = readtable(reader_file);
group_idx = reader_data.ms_b_twre_swe_rawfinal < median(reader_data.ms_b_twre_swe_rawfinal);
bad_readers = reader_data.participant_id(group_idx);
% just to be sure that our group idx matches the order in folder_names 
group_idx = cellfun(@(x) contains(x, bad_readers), folder_names)';

%% Run on the first 3 conditions (all oddball)
% (1) pseudo-fonts with word oddball
% (2) pseudo-fonts with non-word oddball 
% (3) non-words with word oddball

use_bins = 0;          % indices of bins to include in analysis (the values must be present in the bin column of all DFT/RLS exports)
%use_bins = [0, 2, 3, 7, 9];      % would typically do this though indices of bins to include in analysis (the values must be present in the bin column of all DFT/RLS exports)
use_freqs = [1 2];       % indices of harmonics to include in analysis 
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

%% NOW RUN RCA ON 

use_subs = group_idx;
rca_read1 = rcaSweep...
    (folder_names, use_bins, use_freqs, use_conds, use_trials, use_subs, n_reg, n_comp, data_type, comp_channel, force_reload, return_all);
if print_figures, print('-dpsc', fullfile(data_location, 'figures', 'rca_read1.ps'), '-append'), end

use_subs = ~group_idx;
rca_read2 = rcaSweep...
    (folder_names, use_bins, use_freqs, use_conds, use_trials, use_subs, n_reg, n_comp, data_type, comp_channel, force_reload, return_all);
if print_figures, print('-dpsc', fullfile(data_location, 'figures', 'rca_read2.ps'), '-append'), end
                       
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
    rca_read1(f) = aggregateData(rca_read1(f), keep_conditions, errorType, trialError, doNR);
    rca_read2(f) = aggregateData(rca_read2(f), keep_conditions, errorType, trialError, doNR);
    % called like this:
    % aggregateData(rca_struct, keep_conditions, error_type, trial_wise, do_nr)
    % more info in 'help'
end

%% COMPARE GROUPS
close all
harm_idx = [1, 2];
comp_idx = [1, 2];
plot_projected = false;

for c = 1:length(comp_idx)
    subplot(1,2,c); 
    hold on
    for z = 1:2
        if z == 1
            cur_idx = group_idx;
            cur_struct = rca_read1;
            x_vals = 1:2:(length(harm_idx)*2);
        else
            cur_idx = ~group_idx;
            cur_struct = rca_read2;
            x_vals = 2:2:(length(harm_idx)*2);
        end
        % note: because all of the data are included in the output, 
        % we need to recompute the means and errors
        real_data = cell2mat(arrayfun(@(x) ...
            squeeze(cur_struct(x).subjects.real_signal(11, 1, comp_idx(c), cur_idx, 1)), harm_idx, 'uni', false));
        imag_data = cell2mat(arrayfun(@(x) ...
            squeeze(cur_struct(x).subjects.imag_signal(11, 1, comp_idx(c), cur_idx, 1)), harm_idx, 'uni', false));
        vector_mean = sqrt(nanmean(real_data, 1).^2 + nanmean(imag_data, 1).^2);
        combined_data = permute(cat(3, real_data, imag_data), [1,3,2]); % move harmonics into third dimension
        for h = 1:length(harm_idx)
            vector_err(h,:) = fitErrorEllipse(combined_data(:,:,harm_idx(h)));
        end
        
        projected_data = cell2mat(arrayfun(@(x) ...
            squeeze(cur_struct(x).subjects.proj_amp_signal(11, 1, comp_idx(c), cur_idx, 1)), harm_idx, 'uni', false));
        projected_mean = nanmean(projected_data, 1);
        projected_err = nanstd(projected_data, [], 1)./sqrt(length(find(cur_idx)));
        projected_err = repmat(projected_err', 1, 2);
        
        if plot_projected
            cur_data = projected_mean;
            cur_err = projected_err;
        else
            cur_data = vector_mean;
            cur_err = vector_err;
        end
        bar(x_vals, cur_data, .4)
        for h = 1:length(harm_idx)
            errorbar(x_vals(h), cur_data(h), projected_err(h,1), projected_err(h,2), 'k-', 'linewidth', l_width)
        end
    end
    set(gca, gca_opts{:}, 'xtick', x_vals-.5, 'xticklabel', rca_read1(1).settings.freqLabels(harm_idx));
    xlim([.5, (length(harm_idx)*2)+.5])
    ylim([0, 6]);
    hold off
end


