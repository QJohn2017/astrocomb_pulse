function [out_t, out_f, power_t, power_f, phase_t, phase_f, fwhm] = prop_smf(pulse, avg_power, z)
% Propagate pulse through Single Mode Fiber

% Center wavelength, nm
lambda = 1550; 

% Fiber attenuation, dB/km
% https://www.thorlabs.com/drawings/55837584f9716cd3-A6762EB0-B4E7-966F-803861B89D5515FC/SMF-28-J9-SpecSheet.pdf
atten = 0.18;

% Loss coefficient, 1/m
alpha = atten_to_alpha(atten);

% SMF-28 zero dispersion wavelength and slope
% https://www.corning.com/media/worldwide/coc/documents/Fiber/SMF-28%20Ultra.pdf
lambda0 = 1314; % nm
slope = .092; % ps/(nm^2 km)

% Dispersion and dispersion slope at 1550 nm
% http://mathscinotes.com/2012/04/optical-fiber-dispersion-formula-where-did-this-come-from/
D = slope/4 * (lambda - lambda0^4/lambda^3); % ps/(nm km)
d_D = slope/4 * (1 + 3 * lambda0^4/lambda^4); % ps/(nm^2 km)

% Dispersion Beta parameters, [1/m ps/m ps^2/m ps^3/m]
[beta2, beta3] = D_to_beta(D, d_D);
betap = [0, 0, beta2, beta3];

% Nonlinear coefficient of material, 1/(W m)
% https://ieeexplore.ieee.org/document/7764544
gamma = 0.78E-3;

% Propagate through the SMF
[out_t, out_f, power_t, power_f, phase_t, phase_f, fwhm] = ...
    prop_pulse(pulse, alpha, betap, gamma, avg_power, z);

end
