function [out_t, out_f, power_t, power_f, phase_t, phase_f, fwhm] = ...
    prop_edfa(pulse, avg_power)
% Propagate pulse through EDFA

% Center wavelength, m
lambda = 1550E-9;

% Parameters for both stages

% Fiber attenuation, dB/km
% https://www.thorlabs.com/drawings/55837584f9716cd3-A6762EB0-B4E7-966F-803861B89D5515FC/SMF-28-J9-SpecSheet.pdf
atten = 0.18;

% Loss coefficient, 1/m
alpha = atten_to_alpha(atten);

% Nonlinear index, m^2/W
% https://www.rp-photonics.com/nonlinear_index.html
n2 = 2.5E-20;

% Mode field diameter and effective mode area
% https://www.thorlabs.com/newgrouppage9.cfm?objectgroup_id=949
mfd = 10.4E-6; % m
A_eff = pi * (mfd/2)^2; % m^2

% Nonlinear coefficient, 1/(W m)
gamma = 2*pi / lambda * n2 / A_eff;

% First stage parameters (from Pritel email and guesses)
z = 5; % Fiber length, m
D = 0.0; % Dispersion, ps/(nm km)
d_D = 0.0; % Dispersion slope, ps/(nm^2 km)
[beta2, beta3] = D_to_beta(D, d_D);
betap = [0, 0, beta2, beta3]; % [1/m, ps/m, ps^2/m, ps^3/m]

% Propagate first stage
out_t = prop_pulse(pulse, alpha, betap, gamma, avg_power(1), z);

% Second stage parameters (from Pritel email and guesses)
z = 11; % Fiber length, m
D = 17.0; % Dispersion, ps/(nm km)
d_D = 0.0; % Dispersion slope, ps/(nm^2 km)
[beta2, beta3] = D_to_beta(D, d_D);
betap = [0, 0, beta2, beta3]; % [1/m, ps/m, ps^2/m, ps^3/m]

% Propagate second stage
[out_t, out_f, power_t, power_f, phase_t, phase_f, fwhm] = ...
    prop_pulse(out_t(:, end), alpha, betap, gamma, avg_power(2), z);

end
