% MATLAB Code to Plot Membership Functions for Member of radiation.

% Range for the variable (e.g., signal strength, density, etc.)
X = 0.5:0.01:2.5;

% Membership functions for low and high clusters
low = 1 ./ (1 + exp(-5 * (1.5 - X)));
high = 1 - low;

% Sample noisy data points for each cluster
low_cluster = low + 0.05 * randn(size(X));
high_cluster = high + 0.05 * randn(size(X));

% Plotting the data points
figure;
scatter(X, low_cluster, 'r*'); hold on;
scatter(X, high_cluster, 'b*');

% Plotting the membership functions
plot(X, low, 'r-', 'LineWidth', 2);
plot(X, high, 'b-', 'LineWidth', 2);

% Labels and Legends
xlabel('X');
ylabel('Degree of membership');
title('Fuzzy Membership Functions for Low and High Clusters');
legend('1st cluster', '2nd cluster', 'membership of 1st cluster', 'membership of 2nd cluster', 'Location', 'Best');

grid on;
axis([0.5 2.5 0 1]);

% Annotating 'low' and 'high' regions
text(0.55, 0.9, 'low', 'FontWeight', 'bold');
text(2.3, 0.9, 'high', 'FontWeight', 'bold');