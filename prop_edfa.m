function [out_t, out_f, power_t, power_f, phase_t, phase_f, fwhm] = ...
    prop_edfa(pulse, avg_power, lambda0)
% Propagate pulse through EDFA

% Length of fiber (m) in each step. Basically a guess
z = [8, 10];

% Dispersion Beta Params as [0, 0, beta2, beta3] unit [0 0 ps2/m ps3/m] 
betap = [0, 0, 0, 5.0308e-05;
         0, 0, -0.0217, 5.0308e-05];
     
% Loss coefficient, unit m-1. Leading number is loss in dB/km
alpha = 0.18 / 1000 / 4.34;
%alpha = (-35/z)/1000/4.34;

% Avg Power, W
avg_power = [0.2, avg_power];

% n2 of material m2/W
n2 = 2.5E-20;

% Effective Area of mode, m^2
Aeff = pi * (6.5 / 2)^2 * (1E-6)^2;

% Nonlinear coefficient 1/(W m)
gamma = (n2/Aeff) * (2*pi / (lambda0*1E-9));

disp('Propagating through EDFA');
tic;
% Loop twice becuase of 2 stage edfa
for i = 1:2
    % Peak power in J(?)
%     peak_power = (energy*1E-12) / (pulse_fwhm*1E-12);
    
    % Initialize pulse from previous pulse
    if i == 1
%         u_ini = sqrt(peak_power) * pulse;
        u_ini = pulse;
    else
%         u_ini = sqrt(peak_power) * power_t(:, end);
        u_ini = power_t(:, end);
    end
    
    [out_t, out_f, power_t, power_f, phase_t, phase_f, fwhm] = ...
        prop_pulse(u_ini, alpha, betap(i, :), gamma, avg_power(i), z(i));
end
toc;

end
