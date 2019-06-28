function [out_t, out_f, power_t, power_f, phase_t, phase_f, fwhm] = prop_smf(pulse, avg_power, z)
% Propagate pulse through Single Mode Fiber

% Loss coefficient, unit m-1. Leading number is loss in dB/km
alpha = 0.18 / 1000 / 4.34;

% Dispersion Beta Params as [0, 0, beta2, beta3] unit [0 0 ps2/m ps3/m] 
betap = 1.0 * [0, 0, -0.0217, 5.0308e-05];

% Nonlinear coefficient of material taken from https://ieeexplore.ieee.org/document/7764544
gamma = 0.78*1E-3;

disp('Propagating through SMF')
tic;
[out_t, out_f, power_t, power_f, phase_t, phase_f, fwhm] = ...
    prop_pulse(pulse, alpha, betap, gamma, avg_power, z);
toc;

end
