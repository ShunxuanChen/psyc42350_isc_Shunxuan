% Task 3
% calculate_isc.m
% Purpose: Calculate temporal, dynamic, and spatial ISC
% Author: [Shunxuan Chen]
% Date: [02/13/2025]

%% Setup
clear all; close all; clc;

% Load the toy data
load('../data/toy_data.mat');

% Get dimensions
[n_subjects, n_rois, n_timepoints, n_sessions] = size(roi_data);

%% Calculate mean time-series across sessions
mean_data = mean(roi_data, 4); % Average across sessions

%% Helper function for leave-one-out mean
function group_mean = get_loo_mean(data, subject_idx)
    % Remove current subject and calculate mean
    other_subjects = true(size(data, 1), 1);
    other_subjects(subject_idx) = false;
    group_mean = mean(data(other_subjects, :), 1);
end

%% 1. Temporal ISC

% Initialize matrices
pairwise_temporal_isc = zeros(n_subjects, n_subjects, n_rois);
loo_temporal_isc = zeros(n_subjects, n_rois);

% Calculate pairwise temporal ISC
for roi = 1:n_rois
    for i = 1:n_subjects
        for j = 1:n_subjects
            if i ~= j
                timeseries_i = squeeze(mean_data(i, roi, :));
                timeseries_j = squeeze(mean_data(j, roi, :));
                pairwise_temporal_isc(i, j, roi) = corr(timeseries_i, timeseries_j);
            end
        end
    end
end

% Calculate leave-one-out temporal ISC
for roi = 1:n_rois
    for subj = 1:n_subjects
        % Get current subject's time series
        subj_ts = squeeze(mean_data(subj, roi, :));
        
        % Get mean of other subjects
        others_ts = get_loo_mean(squeeze(mean_data(:, roi, :)), subj);
        
        % Calculate correlation
        loo_temporal_isc(subj, roi) = corr(subj_ts, others_ts');
    end
end

%% 2. Dynamic ISC
window_size = 10;
step_size = 10;
n_windows = floor(n_timepoints/window_size);

% Initialize dynamic ISC matrix
loo_dynamic_isc = zeros(n_subjects, n_rois, n_windows);

% Calculate dynamic ISC
for roi = 1:n_rois
    for subj = 1:n_subjects
        for w = 1:n_windows
            % Get time window indices
            window_start = (w-1)*step_size + 1;
            window_end = min(window_start + window_size - 1, n_timepoints);
            
            % Get current subject's window time series
            subj_ts = squeeze(mean_data(subj, roi, window_start:window_end));
            
            % Get mean of other subjects for this window
            others_ts = get_loo_mean(squeeze(mean_data(:, roi, window_start:window_end)), subj);
            
            % Calculate correlation
            loo_dynamic_isc(subj, roi, w) = corr(subj_ts, others_ts');
        end
    end
end

% Plot dynamic ISC for first subject and ROI
figure('Name', 'Dynamic ISC');
plot(1:n_windows, squeeze(loo_dynamic_isc(1, 1, :)), '-o', 'LineWidth', 2);
xlabel('Time Window');
ylabel('ISC Value');
title('Dynamic ISC - Subject 1, ROI 1');
grid on;

% Save plot
saveas(gcf, '../data/dynamic_isc_plot.png');

%% 3. Spatial ISC

% Calculate temporal mean for each ROI
temporal_mean_data = mean(mean_data, 3); % Average over time
loo_spatial_isc = zeros(n_subjects, 1);

% Calculate spatial ISC
for subj = 1:n_subjects
    % Get current subject's spatial pattern
    subj_pattern = temporal_mean_data(subj, :);
    
    % Get mean pattern of other subjects
    others_pattern = get_loo_mean(temporal_mean_data, subj);
    
    % Calculate correlation
    loo_spatial_isc(subj) = corr(subj_pattern', others_pattern');
end

%% Save results
save('../data/isc_results.mat', 'pairwise_temporal_isc', 'loo_temporal_isc', ...
    'loo_dynamic_isc', 'loo_spatial_isc');

%% Display summary statistics
fprintf('\nSummary Statistics:\n');
fprintf('Temporal ISC (LOO) - Mean: %.3f, Std: %.3f\n', ...
    mean(loo_temporal_isc(:)), std(loo_temporal_isc(:)));
fprintf('Dynamic ISC - Mean: %.3f, Std: %.3f\n', ...
    mean(loo_dynamic_isc(:)), std(loo_dynamic_isc(:)));
fprintf('Spatial ISC - Mean: %.3f, Std: %.3f\n', ...
    mean(loo_spatial_isc), std(loo_spatial_isc));

% Additional visualization for verification
figure('Name', 'ISC Results Overview');

% Temporal ISC
subplot(3,1,1);
imagesc(mean(pairwise_temporal_isc, 3));
colorbar;
title('Mean Pairwise Temporal ISC across ROIs');
xlabel('Subject');
ylabel('Subject');

% Dynamic ISC
subplot(3,1,2);
imagesc(squeeze(mean(loo_dynamic_isc, 1)));
colorbar;
title('Dynamic ISC across ROIs and Windows');
xlabel('Time Window');
ylabel('ROI');

% Spatial ISC
subplot(3,1,3);
bar(loo_spatial_isc);
title('Spatial ISC by Subject');
xlabel('Subject');
ylabel('ISC Value');

saveas(gcf, '../data/isc_overview.png');