clc
clear variables
close all

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


%% State matrix representation

global nx nu
nx = 4;
nu = 1;

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
funcs.phi  = @(x)    1/2 * (x - x_final)' * funcs.P * (x - x_final);
funcs.L    = @(x, u) 1/2 * (x'*funcs.Q*x + u'*funcs.R*u);
funcs.Lx   = @(x, u) funcs.Q*x;
funcs.Lu   = @(x, u) funcs.R*u;
funcs.phix = @(x, u) funcs.P * (x - x_final);


%% Iterative Procedure

N = 100;
dt = (t_final - t_initial) / (N-1);

% Initial conditions for the minimization problem
u = consts.T * ones(N-1, 1);
x = zeros(N, nx);
x(1, :) = x_initial';

for ii = 1:N-1
    x(ii+1, :) = x(ii, :) + dt * compute_x_dot(x(ii, :), u(ii, :), consts)';
end

z0 = combine_to_z(x, u);

[z, J] = fmincon( ...
    @(z) compute_J_gradJ(z, dt, consts, funcs), ...
    z0, ...
    [], [], ...
    [], [], ...
    [], combine_to_z(+Inf(N, nx), consts.T * ones(N-1, 1)), ...
    @(z) compute_C_gradC(z, dt, x_initial, consts, funcs), ...
    optimoptions('fmincon', ...
        'MaxFunctionEvaluations', 50000, ...
        'SpecifyObjectiveGradient', true, ...
        'Display', 'iter') ...
    );

[x, u] = extrapolate_from_z(z);


%% Plots

reset(0);
set(0, 'DefaultFigureNumberTitle', 'off');
set(0, 'DefaultFigureWindowStyle', 'docked');
set(0, 'defaultaxesfontsize', 12);
set(0, 'DefaultLineLineWidth', 2);

[figure_trajectory] = plot_trajectory(x', [u; 0]', x_final, t_initial:dt:t_final, 'Direct transcription (direct)');

%%
exportgraphics(figure_trajectory, 'resources/exports/direct_trajectory.jpg', 'Resolution', 600)


%% Functions

function [x_dot] = compute_x_dot(x, u, consts)

x_dot = zeros(length(x), 1);

h = x(1); %#ok<NASGU>
v = x(2);
m = x(3);
gamma = x(4);

x_dot(1) = v * sin(gamma);
x_dot(2) = u/m - consts.Cd*v^2/m - consts.g * sin(gamma);
x_dot(3) = - consts.alpha * u;
x_dot(4) = consts.Cl*v/m - consts.g * cos(gamma) / v;

end


function [J, gradJ] = compute_J_gradJ(z, dt, ~, funcs)

[x, u] = extrapolate_from_z(z);
    
J = 0;
for ii = 1:size(x, 1)-1
    J = J + dt * funcs.L(x(ii, :)', u(ii, :));
end
J = J + funcs.phi(x(end, :)');


% Gradient of the objective function
if (nargout > 1)
    
    gradJx = zeros(size(x));
    gradJu = zeros(size(u));

    for ii = 1:size(x, 1)-1
        gradJx(ii, :) = dt * funcs.Lx(x(ii, :), u(ii, :));
        gradJu(ii, :) = dt * funcs.Lu(x(ii, :), u(ii, :));    
    end

    gradJx(end, :) = funcs.phix(x(end, :)')';

    gradJ = combine_to_z(gradJx, gradJu);

end

end


function [C, C_vector, gradC, gradC_vector] = compute_C_gradC(z, dt, x_initial, consts, funcs)

[x, u] = extrapolate_from_z(z);

[N, nx] = size(x);
[~, nu] = size(u);

% Constraint function
C = [];
C_vector = zeros(N, nx);
C_vector(1, :) = x_initial' - x(1, :);

for ii = 1:N-1
    C_vector(ii+1, :) = x(ii, :) + dt * compute_x_dot(x(ii, :), u(ii, :), consts)' - x(ii + 1, :);
end
C_vector = combine_to_z(C_vector);


if (nargout > 2)
    gradC = [];
    gradC_vector = zeros(nx, N*(nu+nx) + nx);
    gradC_vector(1:nx, 1:nx) = - eye(nx);

    for ii = 1:N-1
        gradC_vector((1 + nx +(ii - 1)*nx):((nx + ii*nx)), (1 +(ii - 1)*(nx+nu)):((ii - 1)*(nx+nu) + nx)) = eye(nx) + dt * funcs.A(x(ii, :), u(ii, :));
        gradC_vector((1 + nx +(ii - 1)*nx):((nx + ii*nx)), (1 +(ii)*(nx+nu)):((ii)*(nx+nu) + nx)) = - eye(nx);
        gradC_vector((1 + nx +(ii - 1)*nx):((nx + ii*nx)), (1 +(ii - 1)*(nx+nu) + nx):((ii - 1)*(nx+nu) + nx + nu)) = + dt * funcs.B(x(ii, :), u(ii, :));
    end

    gradC_vector = gradC_vector.';
    
end

end


function z = combine_to_z(x, u)

global nu

if (nargin == 1)
    z = x';
    z = z(:);
    return
end

z = [x'; [u' zeros(nu, 1)]];
z = z(1:end - nu);

end


function [x, u] = extrapolate_from_z(z)
% Extract states and control from z

global nx nu

z = [z zeros(1, nu)];
z = reshape(z, nx + nu, [])';

x = z(:, 1:nx);
u = z(1:end-1, nx+1:end);

end

