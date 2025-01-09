clc
clear variables
close all

%% Operating point and system matrices

u0 = 5; %[V]
[x_eq, u_eq] = operating_point(u0);
[A, B, C, D] = ABCD(x_eq, u_eq);
G = tf(ss(A, B, C, D));


%% Controller (LQR)

% Bryson’s rule
Q = diag([1 0 1 0] ./ ([0.3 1 2 25].^2));
R = diag(0.5 ./ 1^2);

K_LQR = lqr(A, B, Q, R);


%% Estimator (EKF)

% Larger q -> I trust the measurements more
Q_KF = diag([0 1e-1 1 1 0]);
R_KF = diag(1e-5);
N_KF = 0;

P0 = diag([1e-2 1 1 1 1]);
x0 = zeros(5, 1);
x0(5) = 10;


%% Simulation run & Plots

reset(0);
set(0, 'DefaultFigureNumberTitle', 'off');
set(0, 'DefaultFigureWindowStyle', 'docked');
set(0, 'defaultaxesfontsize', 12);
set(0, 'DefaultLineLineWidth', 2);

figure('Name', 'Simulation results')
tiles = tiledlayout(2, 3, 'TileSpacing', 'tight');

for sim_idx = 1:3

    switch sim_idx
        case 1
        % P0 = diag([1e-2 1 1 1 0.1]);
        % Q_KF = diag([0 1e-1 1 1 0]);
        R_KF = diag(1e-2);
        case 2
        % P0 = diag([1e-2 1 1 1 1]);
        % Q_KF = diag([0 1e-1 1 1 1e-1]);
        R_KF = diag(1e-5);
        case 3
        % P0 = diag([1e-2 1 1 1 10]);
        % Q_KF = diag([0 1e-1 1 1 1]);
        R_KF = diag(1e-8);
    end

    simConfig.StopTime = num2str((x0(5) == 0) * 5 + (x0(5) ~= 0) * 5);
    out = sim("System.slx", simConfig);

    states = out.states;
    covariance = out.covariance;
    noise = out.noise;
    control = out.control;

    titles = {'Position', 'Velocity', 'Current', 'Temperature', 'Parameter (alpha)'};
    ylabels = {'[m]', '[m/s]', '[A]', '[°C]', '[N/A]'};

    colors = ['r', 'g', 'b'];

    for ii = 1:3
        
        a = [1 2 5];
        state_idx = a(ii);
        
        nexttile(tiles, ii);
        hold on
        grid on

        plot(states.time, (states.Data(:, state_idx) - states.Data(:, state_idx + 5)), [colors(sim_idx) '-'], 'DisplayName', ['x' num2str(state_idx) '(t0)=' num2str(x0(state_idx))]);

        title([titles{state_idx} ' error'])
        ylabel(['$x_ ' num2str(state_idx) ' - \hat{x}_ ' num2str(state_idx) ' $ \quad' ylabels{state_idx}], 'Interpreter', 'latex', 'FontSize', 16)
        xlabel('Time [s]')
        ylim('padded')
        legend('Location', 'best')
        
        nexttile(tiles, ii + 3);
        hold on
        grid on
        covariance_tmp = sqrt(abs(covariance.Data(:, sub2ind([5, 5], state_idx, state_idx))));
        plot(covariance.time, +covariance_tmp, [colors(sim_idx) '--'], 'DisplayName', ['R(t0)=' num2str(R_KF)], 'HandleVisibility', 'on');
        plot(covariance.time, -covariance_tmp, [colors(sim_idx) '--'], 'HandleVisibility', 'off');

        title([titles{state_idx} ' covariance'])
        ylabel(['$\sqrt{P_{' num2str(state_idx) num2str(state_idx) '}} $ \quad' ylabels{state_idx}], 'Interpreter', 'latex', 'FontSize', 16)
        xlabel('Time [s]')
        ylim('padded')
        legend('Location', 'best')

    end

end

nexttile(tiles, 1)
hold on
plot(noise.time, noise.Data(:, 1)/2, 'Color', [0 0 0 0.2], 'DisplayName', 'Noise (v)')

%%
% exportgraphics(gcf, 'resources/exports/EKF_results_R.jpg', 'Resolution', 600)