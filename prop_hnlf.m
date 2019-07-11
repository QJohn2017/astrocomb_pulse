function [out_t, out_f, power_t, power_f, phase_t, phase_f, fwhm] = ...
    prop_hnlf(pulse, avg_power, z)
% Propagate pulse through Highly Nonlinear Fiber

% Fiber attenuation, dB/km
% https://fiber-optic-catalog.ofsoptics.com/item/optical--fibers/highly-nonlinear-fiber-optical-fibers1/hnlf-standard-highly-non-linear-fiber-modules
atten = 0.22;

% Loss coefficient, 1/m
alpha = atten_to_alpha(atten);

% Dispersion and dispersion slope of the HNLF
D = -1.8; % ps/(nm km)
d_D = .019; % ps/(nm^2 km)

% Dispersion Beta parameters, [1/m ps/m ps^2/m ps^3/m]
[beta2, beta3] = D_to_beta(D, d_D);
betap = [0, 0, beta2, beta3];

% Nonlinear coefficient, 1/(W m)
gamma = 11.5 * 1E-3;

% Propagate the pulse
[out_t, out_f, power_t, power_f, phase_t, phase_f, fwhm] = ...
    prop_pulse(pulse, alpha, betap, gamma, avg_power, z);

end
