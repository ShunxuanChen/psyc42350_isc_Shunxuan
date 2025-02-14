% Task 5
% calculate_isfc.m
% Purpose: Calculate inter-subject functional connectivity
% Author: [Shunxuan Chen]
% Date: [02/13/2025]

%% Setup
clear all; close all; clc;

% Load the data
load('../data/toy_data.mat');

% Get dimensions
[n_subjects, n_rois, n_timepoints, n_sessions] = size(roi_data);

%% Calculate mean time-series across sessions
mean_data = mean(roi_data, 4); % Average across sessions

%% Calculate leave-one-out ISFC
loo_isfc = zeros(n_subjects, n_rois, n_rois);

for subj = 1:n_subjects
    % Get current subject's time series for all ROIs
    subj_ts = squeeze(mean_data(subj, :, :)); % [n_rois x n_timepoints]
    
    % Get mean time series of all other subjects
    other_subjects = true(n_subjects, 1);
    other_subjects(subj) = false;
    others_mean_ts = squeeze(mean(mean_data(other_subjects, :, :), 1)); % [n_rois x n_timepoints]
    
    % Calculate correlation between all ROI pairs
    % (current subject's ROIs vs mean of other subjects' ROIs)
    loo_isfc(subj, :, :) = corr(subj_ts', others_mean_ts');
end

%% Save results
save('../data/isfc_results.mat', 'loo_isfc');

%% Visualize results
% Plot mean ISFC matrix across subjects
mean_isfc = squeeze(mean(loo_isfc, 1));

figure('Name', 'Mean ISFC Matrix');
imagesc(mean_isfc);
colorbar;
title('Mean Inter-Subject Functional Connectivity');
xlabel('ROI (Other Subjects)');
ylabel('ROI (Target Subject)');
axis square;

% Add ROI labels
xticks(1:n_rois);
yticks(1:n_rois);
xticklabels(arrayfun(@(x) sprintf('ROI %d', x), 1:n_rois, 'UniformOutput', false));
yticklabels(arrayfun(@(x) sprintf('ROI %d', x), 1:n_rois, 'UniformOutput', false));

% Save plot
saveas(gcf, '../data/mean_isfc_plot.png');

%% Additional visualization: ISFC variability
std_isfc = squeeze(std(loo_isfc, 0, 1));

figure('Name', 'ISFC Variability');
imagesc(std_isfc);
colorbar;
title('Standard Deviation of ISFC across Subjects');
xlabel('ROI (Other Subjects)');
ylabel('ROI (Target Subject)');
axis square;

% Add ROI labels
xticks(1:n_rois);
yticks(1:n_rois);
xticklabels(arrayfun(@(x) sprintf('ROI %d', x), 1:n_rois, 'UniformOutput', false));
yticklabels(arrayfun(@(x) sprintf('ROI %d', x), 1:n_rois, 'UniformOutput', false));

saveas(gcf, '../data/std_isfc_plot.png');

%% Display summary statistics
fprintf('\nISFC Summary Statistics:\n');
fprintf('Overall Mean ISFC: %.3f\n', mean(loo_isfc(:)));
fprintf('Overall Std ISFC: %.3f\n', std(loo_isfc(:)));

% Display strongest connections
[max_vals, max_idx] = maxk(mean_isfc(:), 3);
[row_idx, col_idx] = ind2sub([n_rois, n_rois], max_idx);

fprintf('\nTop 3 strongest connections:\n');
for i = 1:3
    fprintf('ROI %d -- ROI %d: %.3f\n', row_idx(i), col_idx(i), max_vals(i));
end

%% Create comprehensive figure with multiple panels
figure('Name', 'ISFC Analysis Overview', 'Position', [100 100 1200 400]);

% Panel 1: Mean ISFC
subplot(1,3,1);
imagesc(mean_isfc);
colorbar;
title('Mean ISFC');
xlabel('ROI (Other Subjects)');
ylabel('ROI (Target Subject)');
axis square;

% Panel 2: Standard deviation of ISFC
subplot(1,3,2);
imagesc(std_isfc);
colorbar;
title('ISFC Variability');
xlabel('ROI (Other Subjects)');
ylabel('ROI (Target Subject)');
axis square;

% Panel 3: Distribution of ISFC values
subplot(1,3,3);
histogram(loo_isfc(:), 30, 'Normalization', 'probability');
title('Distribution of ISFC Values');
xlabel('ISFC Value');
ylabel('Probability');

saveas(gcf, '../data/isfc_overview.png');
