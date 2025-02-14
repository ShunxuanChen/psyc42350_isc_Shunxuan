% Task 4
% calculate_intrasubject_isc.m
% Purpose: Calculate intra-subject correlation between sessions
% Author: [Your Name]
% Date: [Current Date]

%% Setup
clear all; close all; clc;

% Load the data
load('../data/toy_data.mat');

% Get dimensions
[n_subjects, n_rois, n_timepoints, n_sessions] = size(roi_data);

%% Calculate intra-subject ISC
intrasubject_temporal_isc = zeros(n_subjects, n_rois);

for subj = 1:n_subjects
    for roi = 1:n_rois
        % Extract time series for both sessions
        session1_ts = squeeze(roi_data(subj, roi, :, 1));
        session2_ts = squeeze(roi_data(subj, roi, :, 2));
        
        % Calculate correlation between sessions
        intrasubject_temporal_isc(subj, roi) = corr(session1_ts, session2_ts);
    end
end

%% Save results
save('../data/intrasubject_isc_results.mat', 'intrasubject_temporal_isc');

%% Visualize results
figure('Name', 'Intra-subject ISC Results');

% Create heatmap of results
subplot(2,1,1);
imagesc(intrasubject_temporal_isc);
colorbar;
title('Intra-subject ISC by Subject and ROI');
xlabel('ROI');
ylabel('Subject');

% Plot mean ISC across subjects for each ROI
subplot(2,1,2);
mean_isc = mean(intrasubject_temporal_isc, 1);
std_isc = std(intrasubject_temporal_isc, 0, 1);
bar(mean_isc);
hold on;
errorbar(1:n_rois, mean_isc, std_isc, 'k.', 'LineWidth', 1.5);
title('Mean Intra-subject ISC by ROI');
xlabel('ROI');
ylabel('Mean ISC');

% Save plot
saveas(gcf, '../data/intrasubject_isc_plot.png');

%% Display summary statistics
fprintf('\nIntra-subject ISC Summary Statistics:\n');
fprintf('Overall Mean ISC: %.3f\n', mean(intrasubject_temporal_isc(:)));
fprintf('Overall Std ISC: %.3f\n', std(intrasubject_temporal_isc(:)));
fprintf('\nMean ISC by ROI:\n');
for roi = 1:n_rois
    fprintf('ROI %d: %.3f (±%.3f)\n', roi, mean_isc(roi), std_isc(roi));
end