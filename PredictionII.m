% MATLAB Code to Plot Predicted vs Measured Air Temperature

% Time axis (synthetic data for demonstration)
time = 0:1:3000;

% Generate sample data for measured and predicted temperatures
measured_temp = 15 + 10 * sin((time/300)) + 2 * randn(size(time));
predicted_temp = measured_temp + 0.5 * sin((time/200));

% Plotting
figure;
plot(time/300, predicted_temp, 'b-', 'LineWidth', 1.5); hold on;
plot(time/300, measured_temp, 'g-', 'LineWidth', 1.5);

% Labels and Legends
xlabel('Time (s)');
ylabel('Temperature (°C)');
title('Predicted vs Measured Air Temperature');
legend('Predict air-temperature', 'Measured air-temperature', 'Location', 'Best');

grid on;