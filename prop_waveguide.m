function [out_t, out_f, power_t, power_f, phase_t, phase_f, fwhm] = ...
    prop_waveguide(pulse, avg_power, z, lambda0)
% Propagate pulse through Aluminum Nitride Waveguide

% Loss coefficient, unit m-1. Leading number is loss in dB/km
alpha = 0.1 / 1000 / 4.34;

% Dispersion Beta Params as [0, 0, beta2, beta3] unit [0 0 ps2/m ps3/m] 
betap = [0, 0, -4.38E-2, -2.41e-4, 3.16e-6];
%betap = [0, 0, -0.0217, 5.0308e-05];

% n2 of materialm m2/W
n2 = 3E-19;

% Effective Area of mode, m^2
Aeff = 1E-6*1E-6;

% Nonlinearity coefficient 1/(W m)
gamma = 2*pi * n2 / Aeff / lambda0 * 1e9;

% Loss after compressor, etc
fiber_loss = .7;

% Transmission from lens to waveguide
trans = .2;

% Apply losses to average power
avg_power = fiber_loss * sqrt(trans) * avg_power;

disp('Propagating through Waveguide');
tic;
[out_t, out_f, power_t, power_f, phase_t, phase_f, fwhm] = ...
    prop_pulse(pulse, alpha, betap, gamma, avg_power, z);
toc;

end
