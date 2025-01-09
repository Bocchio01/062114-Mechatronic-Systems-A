clc
clear variables
% close all

%% Problem data

consts = struct( ...
    'G', 6.67 * 1e-11, ...
    'g', 9.81, ...
    'M', 5.972 * 1e24, ...
    'R', 6378 * 1e3, ...
    'Cl', 1/2 * 1.29 * 2.05, ...
    'Cd', 1/2 * 1.29 * 0.78, ...
    'T', (13000 * 2 + 1750 * 3) * 1e3, ...
    'alpha', NaN);

% Time settings
t_initial = 0;
t_final = 120;

time_of_u = linspace(t_initial, t_final, 100);

% States settings
x_initial = [0    50   2050e3 deg2rad(89)]';
x_final   = [50e3 1200 880e3  deg2rad(30)]';

consts.alpha = -(x_final(3) - x_initial(3)) / (t_final - t_initial) / (0.85*consts.T);


%% State matrix representation

funcs = struct();

% Q, R, P matrices
funcs.Q = 0;
funcs.R = 0.9 * 1e-6;
funcs.P = 1e+3 * diag([3e-1 1e-0 1e-5 1e-3]);

funcs.A = @(x, u) [
    0 sin(x(4)) 0 x(2)*cos(x(4));
    0 -2*consts.Cd*x(2)/x(3) -(u-consts.Cd*x(2)^2)/x(3)^2 -consts.g*cos(x(4));
    0 0 0 0;
    0 consts.Cl/x(3)+consts.g*cos(x(4))/x(2)^2 -consts.Cl*x(2)/x(3)^2 consts.g*sin(x(4))/x(2)
    ];

funcs.B = @(x, u) [
    0
    1/x(3)
    -consts.alpha
    0
    ];

% Cost function terms
funcs.phi = @(x)    1/2 * (x - x_final)' * funcs.P * (x - x_final);
funcs.L   = @(x, u) 1/2 * (x'*funcs.Q*x + u'*funcs.R*u);
% funcs.L   = @(x, time_of_x, u, time_of_u) 2*...
%     sum(rowfun(@(h, v, m, alpha) 1/2 * [h, v, m, alpha] * funcs.Q * [h, v, m, alpha]', array2table(x(1:end-1, :)), "OutputFormat", "uniform")' * diff(time_of_x)) + ...
%     sum(rowfun(@(T) 1/2 * (T) * funcs.R * (T)', array2table(u(1:end-1, :)), "OutputFormat", "uniform")' * diff(time_of_u)');


%% Iterative Procedure

u = consts.T * ones(length(time_of_u), 1);

step = consts.T * 1e-5;
eps = 1e-2;

J = NaN;

for iter_idx = 1:1000

    % Integrate EoM considering the control resulted from previous cycle
    [time_of_x, x] = ode45(@(t, x) ...
        compute_x_dot(t, x, time_of_u, u, consts), ...
        [t_initial t_final], ...
        x_initial, ...
        odeset('RelTol', 1e-4, 'AbsTol', 1e-6));

    % Compute the adjoint vectors trajectory (the matrix of lambdas)
    x_time_final = x(end, :)';
    lambda_final = (funcs.P * (x_time_final - x_final));

    [time_of_lambdas, lambdas] = ode45(@(t, lambdas) ...
        compute_lambdas_dot(t, lambdas, time_of_x, x, time_of_u, u, funcs), ...
        [t_final t_initial], ...
        lambda_final, ...
        odeset('RelTol', 1e-4, 'AbsTol', 1e-8));

    % We need to align the ode45 results for lambda on the same time grid
    % used by the states x
    lambdas = interp1(time_of_lambdas, lambdas, time_of_x);
    time_of_lambdas = time_of_x;

    % Compute Hamiltonian energetic sensitivity to control input dHdu
    dHdu = compute_dHdu(time_of_x, x, time_of_u, u, time_of_lambdas, lambdas, funcs);

    % Compute the cost function J
    L = 0;
    for t_idx = 1:length(time_of_x) - 1
        L = L + x(t_idx, :)*funcs.Q*x(t_idx, :)' * diff(time_of_x(t_idx:t_idx+1));
    end
    for t_idx = 1:length(time_of_u) - 1
        L = L + u(t_idx, :)*funcs.R*u(t_idx, :)' * diff(time_of_u(t_idx:t_idx+1));
    end

    J(iter_idx) = funcs.phi(x_time_final) + L;

    % Exit condition
    if norm(dHdu) < eps
        break;
    end

    % Adjust control for next iteration
    u = u - step * interp1(time_of_x, dHdu, time_of_u)';
    u(u > consts.T) = consts.T;

    if (rem(iter_idx, 20) == 0)
        clc
        disp(['Iteration #' num2str(iter_idx)])
    end

end

time = time_of_u;
x = interp1(time_of_x, x, time_of_u);
u = interp1(time_of_u, u, time_of_u);


%% Plots

reset(0);
set(0, 'DefaultFigureNumberTitle', 'off');
set(0, 'DefaultFigureWindowStyle', 'docked');
set(0, 'defaultaxesfontsize', 12);
set(0, 'DefaultLineLineWidth', 2);

% [figure_optimization] = plot_optimization(x_TMP, u_TMP, J_TMP, time, true);
[figure_trajectory] = plot_trajectory(x', u, x_final, time, 'Gradient descend method (indirect)');

%%
% exportgraphics(figure_trajectory, 'resources/exports/indirect_trajectory.jpg', 'Resolution', 600)


%% Functions

function [x_dot] = compute_x_dot(t, x, time_of_u, u, consts)

x_dot = zeros(length(x), 1);

h = x(1); %#ok<NASGU>
v = x(2);
m = x(3);
gamma = x(4);

u = interp1(time_of_u, u, t);

x_dot(1) = v * sin(gamma);
x_dot(2) = u/m - consts.Cd*v^2/m - consts.g * sin(gamma);
x_dot(3) = - consts.alpha * u;
x_dot(4) = consts.Cl*v/m - consts.g * cos(gamma) / v;


end


function lambdas_dot = compute_lambdas_dot(t, lambdas, time_of_x, x, time_of_u, u, funcs)

Lx  = @(x, u) funcs.Q*x;

x = interp1(time_of_x, x, t)';
u = interp1(time_of_u, u, t)';

lambdas(isnan(lambdas)) = 0;
lambdas(isinf(lambdas)) = 1e20;

lambdas_dot = -(Lx(x, u) + (lambdas' * funcs.A(x, u))');

if any(~isfinite(lambdas_dot))
    lambdas_dot = 1e20 .* sign(lambdas_dot) ;
end

end


function dHdu = compute_dHdu(time_of_x, x, time_of_u, u, time_of_lambdas, lambdas, funcs)

Lu  = @(x, u) funcs.R*u;

dHdu = zeros(length(lambdas), 1);
x = interp1(time_of_x, x, time_of_lambdas);
u = interp1(time_of_u, u, time_of_lambdas);

for time_idx = 1:length(time_of_lambdas)    
    dHdu(time_idx) = Lu(x(time_idx, :), u(time_idx, :)) + lambdas(time_idx, :) * funcs.B(x(time_idx, :), u(time_idx, :));
end

end

