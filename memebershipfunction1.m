% MATLAB Code to Plot Membership Functions for Clusters For Internal temperature

% Temperature range
T = 18:0.1:25;

% Membership functions for hot and cold clusters
hot = 1 ./ (1 + exp(-5 * (21 - T)));
cold = 1 - hot;

% Sample noisy data points for each cluster
hot_cluster = hot + 0.05 * randn(size(T));
cold_cluster = cold + 0.05 * randn(size(T));

% Plotting the data points
figure;
scatter(T, hot_cluster, 'r*'); hold on;
scatter(T, cold_cluster, 'b*');

% Plotting the membership functions
plot(T, hot, 'r-', 'LineWidth', 2);
plot(T, cold, 'b-', 'LineWidth', 2);

% Labels and Legends
xlabel('Temperature (°C)');
ylabel('Degree of membership');
title('Fuzzy Membership Functions for Hot and Cold Clusters');
legend('1st cluster', '2nd cluster', 'membership of 1st cluster', 'membership of 2nd cluster', 'Location', 'Best');

grid on;
axis([18 25 0 1]);

% Annotating 'hot' and 'cold' regions
text(18.5, 0.9, 'hot', 'FontWeight', 'bold');
text(24, 0.9, 'cold', 'FontWeight', 'bold');
