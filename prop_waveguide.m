function [out_t, out_f, power_t, power_f, phase_t, phase_f, fwhm] = ...
    prop_waveguide(pulse, avg_power, z)
% Propagate pulse through Aluminum Nitride Waveguide

% Center wavelength, nm
lambda = 1550E-9;

% Waveguide attenuation, dB/km
atten = 0.1;

% Loss coefficient, m-1
alpha = atten_to_alpha(atten);

% Dispersion Beta Params as [0, 0, beta2, beta3] unit [0 0 ps2/m ps3/m] 
betap = [0, 0, -4.38E-2, -2.41e-4, 3.16e-6];

% n2 of material m^2/W
n2 = 3E-19;

% Effective Area of mode, m^2
Aeff = 1E-6*1E-6;

% Nonlinearity coefficient 1/(W m)
gamma = 2*pi / lambda * n2 / Aeff;

disp('Propagating through Waveguide');
tic;
[out_t, out_f, power_t, power_f, phase_t, phase_f, fwhm] = ...
    prop_pulse(pulse, alpha, betap, gamma, avg_power, z);
toc;

end
