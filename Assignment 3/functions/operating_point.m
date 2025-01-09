function [x_eq, u_eq] = operating_point(u_star)

% Load system parameters
load("parameters.mat"); %#ok<LOAD>

x_eq = zeros(4, 1);
u_eq = zeros(1, 1);

x_eq(3) = u_star / R;
x_eq(2) = 0;

root_x1 = roots([k3 0 k1 -alpha0*x_eq(3)]);
root_x1 = root_x1((imag(root_x1) == 0) & (real(root_x1) > 0));
assert(~isempty(root_x1), 'No real solution to piston displacement');
x_eq(1) = root_x1(1);
x_eq(4) = (R * x_eq(3)^2) / h + T_env;

u_eq(1) = u_star;


end

