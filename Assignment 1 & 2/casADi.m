clc
clear variables
close all

import casadi.*

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

% States settings
x_initial = [0    50   2050e3 deg2rad(89)]';
x_final   = [50e3 1200 880e3  deg2rad(30)]';

consts.alpha = -(x_final(3) - x_initial(3)) / (t_final - t_initial) / (0.85*consts.T);


%% Optimization parameters

N = 100;
time = linspace(t_initial, t_final, N+1);
dt = diff(time(1:2));

% Cost function weights
Q = 0; % [1e-4 1e-2 1e-6 1-10]';
R = 0.9*1e-6;
P = [1e-2 1e-0 1e-3 1e-7]';

% Q = 0;
% R = 0.9 * 1e-6;
% P = 1e+3 * diag([1e-1 1e-0 1e-5 1e-3]);

x_desired = repmat(x_final, [1, N+1]);

%% Model of the system

h = SX.sym('h');
v = SX.sym('v');
m = SX.sym('m');
gamma = SX.sym('gamma');

T = SX.sym('T');

state = [h; v; m; gamma];
u = T;

nx = length(state);
nu = length(u);

G = [
    v * sin(gamma); ...
    T / m - consts.Cd*v^2 / m - consts.g * sin(gamma); ...
    -consts.alpha * T; ...
    consts.Cl*v / m - consts.g * cos(gamma) / v
    ];

f = Function('f', {state, u}, {G});


%% Optimization stuff

% Optimization variables
opti = Opti();
U = opti.variable(nu, N);
X = opti.variable(nx, N+1);

% Create the cost function
obj = 0;

X_guess = zeros(4, N+1);
for x_idx = 1:nx
    X_guess(x_idx, :) = linspace(x_initial(x_idx), x_final(x_idx), N+1);
end

opt = load("casADi (optimal solution).mat");
opti.set_initial(X, interp1(linspace(0, max(time), size(opt.x_opt, 2)), opt.x_opt', time)');
opti.set_initial(U, interp1(linspace(0, max(time), size(opt.u_opt, 2)), opt.u_opt', time(1:end-1))');

opti.subject_to(X(:, 1) == x_initial)
% opti.subject_to(X(1, :) >= x_initial(1))
% opti.subject_to(X(2, :) >= 1e-3)
% opti.subject_to(X(3, :) >= x_final(3))
% opti.subject_to(X(3, :) <= x_initial(3))
% opti.subject_to(X(4, :) <= x_initial(4))
% opti.subject_to(X(4, :) >= x_final(4))
opti.subject_to(U(1, :) <= consts.T)

for i=1:N
    u_i = U(:, i);
    x_i = X(:, i);
    obj = obj + (sumsqr(x_i .* Q)  + sumsqr(u_i .* R)) * dt;

    % Continuity constraint
    x_dot = f(x_i, u_i);
    opti.subject_to(X(:, i+1) == X(:, i) + x_dot * dt)
end

% Final position constraint
obj = obj + sumsqr((X(:, N) - x_final) .* P);


opts = struct;
% opts.ipopt.print_level = 1;
opts.print_time = 0;

opti.solver('ipopt', opts);
opti.minimize(obj);
opti.solve;

x_opt = opti.value(X);
u_opt = opti.value(U);

save('casADi (optimal solution)', 'time', 'x_opt', 'u_opt')


%% Plots

reset(0);
set(0, 'DefaultFigureNumberTitle', 'off');
% set(0, 'DefaultFigureWindowStyle', 'docked');
set(0, 'defaultaxesfontsize', 15);
set(0, 'DefaultLineLineWidth', 1);

[figure_trajectory] = plot_trajectory(x_opt, [u_opt u_opt(end)], x_final, time, 'casADi solver (direct)');
exportgraphics(figure_trajectory, 'submitted/casADi.jpg', 'Resolution', 600)