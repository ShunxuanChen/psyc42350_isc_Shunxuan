% Task 2
% generate_toy_data.m
% Purpose: Generate simulated fMRI and behavioral data
% Author: [Your Name]
% Date: [Current Date]

%% Clear workspace and set random seed for reproducibility
clear all; close all; clc;
rng(42); % Set random seed

%% Set parameters
n_subjects = 20;    % Number of subjects
n_rois = 5;        % Number of ROIs
n_timepoints = 50;  % Number of time points
n_sessions = 2;     % Number of scan sessions per subject

%% Generate ROI data
% Initialize 4-D matrix (subjects x ROIs x timepoints x sessions)
roi_data = zeros(n_subjects, n_rois, n_timepoints, n_sessions);

% Fill with random data following a realistic pattern
for subj = 1:n_subjects
    for sess = 1:n_sessions
        for roi = 1:n_rois
            % Generate time series with some temporal autocorrelation
            baseline = randn(1) * 0.5; % Random baseline for each ROI
            signal = smooth(randn(n_timepoints, 1)) * 0.5 + baseline;
            roi_data(subj, roi, :, sess) = signal;
        end
    end
end

%% Generate behavioral scores
% Generate behavior scores that follow a normal distribution
behavior = randn(1, n_subjects);

%% Save the data
save('../data/toy_data.mat', 'roi_data', 'behavior');

%% Visualize the data to verify
% Plot ROI time series for first subject, first session
figure('Name', 'Example ROI Time Series');
for roi = 1:n_rois
    subplot(n_rois, 1, roi);
    plot(squeeze(roi_data(1, roi, :, 1)));
    title(['ROI ' num2str(roi) ' - Subject 1, Session 1']);
    xlabel('Time Points');
    ylabel('Signal');
end

% Plot behavioral scores distribution
figure('Name', 'Behavioral Scores Distribution');
histogram(behavior, 10);
title('Distribution of Behavioral Scores');
xlabel('Score');
ylabel('Frequency');

fprintf('Data generation complete!\n');
fprintf('Generated roi_data with dimensions: [%d, %d, %d, %d]\n', ...
    size(roi_data, 1), size(roi_data, 2), size(roi_data, 3), size(roi_data, 4));
fprintf('Generated behavior vector with length: %d\n', length(behavior));