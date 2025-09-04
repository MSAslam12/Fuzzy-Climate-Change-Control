% MATLAB Code to Plot Membership Functions for Internal humidity.

% Humidity index range
Hi = 50:0.1:75;

% Membership functions for dry and wet clusters
dry = 1 ./ (1 + exp(-5 * (60 - Hi)));
wet = 1 - dry;

% Sample noisy data points for each cluster
dry_cluster = dry + 0.05 * randn(size(Hi));
wet_cluster = wet + 0.05 * randn(size(Hi));

% Plotting the data points
figure;
scatter(Hi, dry_cluster, 'r*'); hold on;
scatter(Hi, wet_cluster, 'b*');

% Plotting the membership functions
plot(Hi, dry, 'r-', 'LineWidth', 2);
plot(Hi, wet, 'b-', 'LineWidth', 2);

% Labels and Legends
xlabel('Hi');
ylabel('Degree of membership');
title('Fuzzy Membership Functions for Dry and Wet Clusters');
legend('1st cluster', '2nd cluster', 'membership of 1st cluster', 'membership of 2nd cluster', 'Location', 'Best');

grid on;
axis([50 75 0 1]);

% Annotating 'dry' and 'wet' regions
text(50.5, 0.9, 'dry', 'FontWeight', 'bold');
text(73, 0.9, 'wet', 'FontWeight', 'bold');
