% MATLAB Code to Plot Temperature Variations Over Time

% Time axis (in minutes)
time = 0:1:750;

% Sample data for temperatures
outdoor_temp = 10 + 15 * sin(time / 120) + 5 * randn(size(time)); % Outdoor temperature (blue)
indoor_temp = 10 + 5 * sin(time / 200) + 2 * randn(size(time));   % Indoor temperature (black)

% Plotting
figure;
plot(time/70, outdoor_temp, 'b-', 'LineWidth', 1.5); hold on;
plot(time/70, indoor_temp, 'k-', 'LineWidth', 1.5);

% Labels and Title
xlabel('Time (minutes)');
ylabel('Temperature (°C)');
title('Temperature Variations Over Time');

grid on;
axis([0 10 -25 35]);
