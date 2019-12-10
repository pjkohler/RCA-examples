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
    code_folder = '/Users/lindseyhasak/code';
    addpath(genpath(sprintf('%s/git/export_fig',code_folder)),'-end');
    addpath(genpath(sprintf('%s/git/rcaBase',code_folder)),'-end');
    addpath(genpath(sprintf('%s/git/mrC',code_folder)),'-end');
    addpath(genpath(sprintf('%s/git/schlegel/matlab_lib/misc',code_folder)),'-end');
    addpath(genpath(sprintf('%s/git/schlegel/matlab_lib/figure',code_folder)),'-end');
else
end


%% Set up some analysis parameters and variables

saveData      = true;                                                                             % set to true if you want to automatically save the output of rcaSweep into the 'outputs' directory
printFigures  = true;                                                                             % set to true if you want to automatically print the figures and save them into the'outputs' directory
save_name  = 'rcaData_Lochy_Cond3';    % name of output                                                         % set to the desired ouput file name - this will contain one big data struct with the results from rcaSweep

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

binsToUse     = 1:10;          %Sweep             % indices of bins to include in analysis (the values must be present in the bin column of all DFT/RLS exports)
%binsToUse     = 0;            %Would typically do this though           % indices of bins to include in analysis (the values must be present in the bin column of all DFT/RLS exports)
freqsToUse    = [1 2];                   % indices of frequencies to include in analysis (the values must be present in the frequency column of all DFT/RLS exports)... these freq were set when importing eegssn
% chose freq 1 and 2, because 3Hz is carrier, so you want nf1 clean
condsToUse    = [1];                      % if you want to include all conditions, create a vector here listing all condition numbers
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
    
else % if it does exist, print warning and give options for overwriting / re-running with different filename
    fprintf('\n\nWARNING: The file %s.mat already exists for this data directory!\n\n',save_path);
    
    if getYN('Would you like to run rcaSweep anyway and OVERWRITE this file?') % do you want to overwrite and run again?
        launchAnalysis = true;
        
    else % if you don't want to do that, would you like to use a new filename instead?
        if getYN('Would you like to run rcaSweep and save the output using a NEW filename?')
            newFileName = input(sprintf('Please type in a new filename and press enter (file already on disk = %s):  ',save_name),'s');
            save_path = [data_location,'/',newFileName];
            launchAnalysis = true;
        end
    end
end

% Now we can do the actual analysis...
if launchAnalysis
    
    rcaStruct = rcaSweep(folder_names,binsToUse,freqsToUse,condsToUse,[],nReg,nComp,dataType,chanToCompare,[],rcPlotStyle,forceSourceData);
    % rcaSweep takes inputs like this:
    % rcaSweep(PATHNAMES,[BINSTOUSE],[FREQSTOUSE],[CONDSTOUSE],[TRIALSTOUSE],[NREG],[NCOMP],[DATATYPE],[COMPARECHAN],[SHOW],[RCPLOTSTYLE])
    % help rcaSweep has more information...
    
    if printFigures, print('-dpsc',[save_path,'.ps'],'-append'), end
    if saveData, save(save_path,'rcaStruct'), end

else % if you in fact *don't* want to run RCA (you ansered n to all of the above checks), just load the data you already have.
    load(save_path);
end

% THIS SECTION OF THE CODE: 
% returns the structure 'rcaStruct', whose dimensions are like this: 
% There will be as many 'structs' as HARMONICS you have analysed in the
% input data
% The output in W and A will be nChanels (128) x n RCA components

%% Average RCA data across participants
% Now that we have run our analyses on individual subjects, we can generate
% group-level statistics by looping across these outputs using
% aggregateData.

keepConditions = true;       % if false, will average across conditions
errorType      = 'SEM';      % can also be 'none' or '95CI' for 95% confidence intervals
trialError     = false;      % if true, will compute trial-wise errors and statistics
doNR           = [];         % Naka-Rushton fitting - turned off for now, but takes a logical array structured like harmonic x rc x condition

% Loop over harmonics and aggregate data!
for f = 1:length(freqsToUse)
    
    fprintf('\n ... Harmonic no. %0.0f ...\n',f);
    rcaStruct(f) = aggregateData(rcaStruct(f),keepConditions,errorType,trialError,doNR);
    % Called like this:
    % aggregateData(rca_struct, keep_conditions, error_type, trial_wise, do_nr)
    % more info in 'help'
 
end

% THIS SECTION OF THE CODE: 
% returns the structure 'rcaStructAgg', whose dimensions are like this:
% There will be as many 'structs' as HARMONICS you have analysed in the
% input data
% The outputs in 'mean' are stored like this: 
% bin (plus last one for vector average) x harmonic (always 1) x RCA component (plus last one - comparison electrode) x condition 
% Have a look in the help for more info.

% Save the data - call it what you like
saveFileNamePath = strcat(out_location,'/AggData_Lochy_Cond3');
save(saveFileNamePath,'rcaStruct')

fprintf('\nDone! \n');

%% FIGURE PARAMS
l_width = 2


%% PLOT ERROR ELLIPSE
close all;
figure
% Specify which bin to plot (usually the average)
bin_to_plot = 'ave';
% Specify which RC component to plot
rc_idx = 2;
% Specify which harmonic to plot (always 1)
h_idx = 1;

% Link the bin to plot to the bin index (11 corresponds to the average of
% the 10 useable bins) 
if bin_to_plot == 'ave'
    b_idx = 11;
else
    b_idx = bin_to_plot;
end

% Concatenate the real signal and the imaginary signal together (bin x
% harmonic x subject x condition)
xy_vals = cat(2, squeeze(rcaStruct(h_idx).subjects.real_signal(b_idx,:,rc_idx,:)), ...
                 squeeze(rcaStruct(h_idx).subjects.imag_signal(b_idx,:,rc_idx,:)));
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

%% Pull out phase significance - another way 

% Shows average t_sig for each component for the first harmonic
squeeze(rcaStruct(1).stats.t_sig(end, 1, :));

% Shows average t2_sig for each component for the first harmonic
squeeze(rcaStruct(1).stats.t2_sig(end, 1, :));

% Can also look at this bin-by-bin
squeeze(rcaStruct(1).stats.t2_sig(:, 1, :))


%% FLIP and REORDER RCA TOPOGRAPHIES
h_idx = 1;
rca_new = flipSwapRCA(rcaStruct(h_idx), [2,1], 2);
figure;
subplot(2,2,1); hold on;
mrC.plotOnEgi(rcaStruct(h_idx).W(:,1))
subplot(2,2,2); 
mrC.plotOnEgi(rcaStruct(h_idx).W(:,2))
subplot(2,2,3); 
mrC.plotOnEgi(rca_new.W(:,1))
subplot(2,2,4); 
mrC.plotOnEgi(rca_new.W(:,2))

%% LOAD IN AXX DATA AND PLOT THEM FOR AN FREQ DOMAIN-DEFINED RC
close all;
rc_idx = 2;
full_w = rca_new.W(:,rc_idx);
axx_rca_full = axxRCAmake(folder_names, full_w);

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


for h  = 1:length(freqsToUse)
    for c = 1:length(nComp)
    % Find the mean amplitude of each component 
    curr_amp = rcaStruct(h).mean.amp_signal(end,:,c); %currently just for 
    curr_errlo = rcaStruct(h).stats.amp_lo_err(end,:, c);
    curr_errhi = rcaStruct(h).stats.amp_up_err(end,:, c);
    end
    plot(comp_n, curr_amp);
end
    
  
 

