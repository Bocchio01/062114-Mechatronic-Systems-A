function [figure_results] = plot_trajectory(x, u, x_final, time, title_string)

figure_results = figure('Name', title_string);
tiles = tiledlayout(5, 2, 'TileSpacing', 'tight');
title(tiles, 'Space Shuttle Simulation');
subtitle(tiles, title_string);

x_position = cumtrapz(time, x(2, :) .* cos(x(4, :)));
y_position = x(1, :);
gamma_position = x(4, :);


%% Tiles setup

% Altitude
t1 = nexttile(tiles, 1);
hold on
grid on

plot(time, x(1, :), '-');
plot(time(end), x_final(1), 'o', 'LineWidth', 1);

title('States and control')
ylim('padded')
ylabel('$h [m]$', 'Interpreter', 'latex', 'FontSize', 14) 
xticklabels("");
set(gca, 'XTick', linspace(0, time(end), 7));

% Velocity
t2 = nexttile(tiles, 3);
hold on
grid on

plot(time, x(2, :), '-');
plot(time(end), x_final(2), 'o', 'LineWidth', 1);

ylim('padded')
ylabel('$v [m/s]$', 'Interpreter', 'latex', 'FontSize', 14) 
xticklabels("");
set(gca, 'XTick', linspace(0, time(end), 7));

% Mass
t3 = nexttile(tiles, 5);
hold on
grid on

plot(time, x(3, :), '-');
plot(time(end), x_final(3), 'o', 'LineWidth', 1);

ylim('padded')
ylabel('$m [kg]$', 'Interpreter', 'latex', 'FontSize', 14) 
xticklabels("");
set(gca, 'XTick', linspace(0, time(end), 7));

% Flight angle
t4 = nexttile(tiles, 7);
hold on
grid on

plot(time, rad2deg(x(4, :)), '-');
plot(time(end), rad2deg(x_final(4)), 'o', 'LineWidth', 1);

ylim('padded')
ylabel('$\gamma [deg]$', 'Interpreter', 'latex', 'FontSize', 14)
xticklabels("");
set(gca, 'XTick', linspace(0, time(end), 7));

% Thrust force
t5 = nexttile(tiles, 9);
hold on
grid on

yline(31250000, '--k', 'Label', 'Maximum thrust')
plot(time, u(1, :), '-');

ylim('padded')
ylabel('$T [N]$', 'Interpreter', 'latex', 'FontSize', 14) 
set(gca, 'XTick', linspace(0, time(end), 7));


linkaxes([t1 t2 t3 t4 t5], 'x')
xlabel('$Time [s]$', 'Interpreter', 'latex', 'FontSize', 14) 
xlim('tight')


nexttile(tiles, 2, [5, 1]);
hold on
grid on

yline(x_final(1), '--k', 'Label', 'Target altitude')
plot(x_position, y_position, '-k')

for ii = floor(linspace(1, length(time), 20))
    quiver(x_position(ii), y_position(ii), cos(gamma_position(ii)), sin(gamma_position(ii)), 15*x(2, ii), 'r');
end

xlim([min(x_position) max(x_position)*1.5])
ylim([min(y_position) max(y_position)*1.5])
title('Trajectory')
xlabel('$x [m]$', 'Interpreter', 'latex', 'FontSize', 14) 
ylabel('$y [m]$', 'Interpreter', 'latex', 'FontSize', 14) 
legend('', 'Trajectory', 'Velocity (scaled 15x)')

end

