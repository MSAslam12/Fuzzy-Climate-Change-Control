% MATLAB Code to Plot Prediction Error in Air Temperature

% Time axis (synthetic data for demonstration)
time = 0:1:3000;

% Generate sample data for measured and predicted temperatures
measured_temp = 14 + 10 * sin((time/300)) + 2 * randn(size(time));
predicted_temp = measured_temp + 0.5 * sin((time/200));

% Calculate error (percentage error)
error_percent = ((predicted_temp - measured_temp) ./ measured_temp) * 100;

% Plotting
figure;
plot(time/300, error_percent, 'b-', 'LineWidth', 1.5);

% Labels and Title
xlabel('Time (s)');
ylabel('Error (%)');
title('Prediction Error in Air Temperature');

grid on;
axis([0 10 -100 125]);
