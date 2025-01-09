clc
clear variables
% close all

%%

C_T = 449; %[J/(kg*°k)] Thermal capacity
h = 80e+3; %[W/k] Thermal dissipation coefficient
T_env = 0; %[] Ambient temperature
m = 0.5; %[kg] Mass of the piston
L0 = 1e-1; %[H] Inductance at a reference temperature
beta1 = 1; %[H/°C] Coefficients for the T-dependent inductance
beta2 = 1; %[H/°C^2] Coefficients for the T-dependent inductance

c = 1e-2; %[] Damping
k1 = 10; %[N/m] Coefficients for the x-dependent stiffness
k3 = 2; %[N/m^3] Coefficients for the x-dependent stiffness

alpha0 = 1.2; %[N/A] Force constant

R = 5; %[ohm] Resistance


%% Updating the .mat file

save('identification/parameters.mat')