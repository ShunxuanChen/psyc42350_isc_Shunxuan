% Task 6
% relate_isc_to_behavior.m
% Purpose: Analyze relationship between ISC and behavioral similarity
% Author: [Shunxuan Chen]
% Date: [02/13/2025]

%% Setup
clear all; close all; clc;

% Load ISC and behavioral data
load('../data/isc_results.mat', 'pairwise_temporal_isc');
load('../data/toy_data.mat', 'behavior');

[n_subjects, ~, n_rois] = size(pairwise_temporal_isc);

%% Calculate behavioral similarity matrix
% Hypothesis: Subjects with more similar behavioral scores (smaller absolute 
% differences in scores) will show more similar brain responses (higher ISC).
% We calculate similarity as the negative absolute difference between scores,
% so higher values indicate more similar behavior.
behavioral_similarity = zeros(n_subjects, n_subjects);
for i = 1:n_subjects
   for j = 1:n_subjects
       behavioral_similarity(i,j) = -abs(behavior(i) - behavior(j));
   end
end

%% Relate ISC to behavioral similarity for each ROI
% Initialize results storage
correlation_results = zeros(n_rois, 1);
p_values = zeros(n_rois, 1);

% For visualization
figure('Name', 'ISC vs Behavioral Similarity');

for roi = 1:n_rois
   % Get ISC values for this ROI
   roi_isc = squeeze(pairwise_temporal_isc(:,:,roi));
   
   % Get upper triangle indices (excluding diagonal)
   upper_idx = triu(true(n_subjects), 1);
   
   % Extract values for correlation
   isc_values = roi_isc(upper_idx);
   behav_values = behavioral_similarity(upper_idx);
   
   % Calculate correlation
   [r, p] = corr(isc_values, behav_values);
   correlation_results(roi) = r;
   p_values(roi) = p;
   
   % Plot relationship for each ROI
   subplot(2,3,roi);
   scatter(behav_values, isc_values, 'filled', 'MarkerFaceAlpha', 0.6);
   hold on;
   
   % Add trend line
   coef = polyfit(behav_values, isc_values, 1);
   x_trend = linspace(min(behav_values), max(behav_values), 100);
   y_trend = polyval(coef, x_trend);
   plot(x_trend, y_trend, 'r--', 'LineWidth', 2);
   
   title(sprintf('ROI %d (r = %.3f, p = %.3f)', roi, r, p));
   xlabel('Behavioral Similarity');
   ylabel('ISC');
   grid on;
end

% Adjust subplot layout
sgtitle('Relationship between ISC and Behavioral Similarity by ROI');

% Save plot
saveas(gcf, '../data/isc_behavior_relationship.png');

%% Display results summary
fprintf('\nCorrelation between ISC and Behavioral Similarity:\n');
fprintf('================================================\n');
for roi = 1:n_rois
   fprintf('ROI %d: r = %.3f (p = %.3f)', roi, correlation_results(roi), p_values(roi));
   if p_values(roi) < 0.05
       fprintf(' *');
   end
   fprintf('\n');
end

%% Additional visualization: Similarity matrices comparison
figure('Name', 'Similarity Matrices Comparison');

% Plot behavioral similarity matrix
subplot(1,2,1);
imagesc(behavioral_similarity);
colorbar;
title('Behavioral Similarity');
xlabel('Subject');
ylabel('Subject');
axis square;

% Plot mean ISC across ROIs
mean_isc = mean(pairwise_temporal_isc, 3);
subplot(1,2,2);
imagesc(mean_isc);
colorbar;
title('Mean ISC across ROIs');
xlabel('Subject');
ylabel('Subject');
axis square;

saveas(gcf, '../data/similarity_matrices.png');

%% Save results
results = struct();
results.behavioral_similarity = behavioral_similarity;
results.correlation_results = correlation_results;
results.p_values = p_values;
save('../data/behavior_isc_relationship.mat', 'results');