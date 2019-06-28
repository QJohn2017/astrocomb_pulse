function [out_t, out_f, power_t, power_f, phase_t, phase_f, fwhm] = prop_hnlf(pulse, avg_power, z)
% Propagate pulse through Highly Nonlinear Fiber

% Loss coefficient, unit m-1. Leading number is loss in dB/km
alpha = 0.22 / 1000 / 4.34;

% Dispersion Beta Params as [0, 0, beta2, beta3] unit [0 0 ps2/m ps3/m] 
betap = [0, 0, 0.0023, -7.2317e-06];

% Nonlinear coefficient 1/(W m)
gamma = 10.6*1E-3;

disp('Propagating through HNLF');
tic;
[out_t, out_f, power_t, power_f, phase_t, phase_f, fwhm] = ...
    prop_pulse(pulse, alpha, betap, gamma, avg_power, z);
toc;

end
