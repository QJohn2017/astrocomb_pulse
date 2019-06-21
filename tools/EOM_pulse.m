function u = EOM_pulse(t, freq_m, P_IM, P_PM)
% Produce a pulse from electro-optic intensity and phase modulators
% Equation taken from https://www.osapublishing.org/ol/abstract.cfm?uri=ol-35-19-3234
%
% INPUT
%
% t - time series
% freq_m - drive frequency in THz
% P_IM - drive power of intensity modulators
% P_PM - drive power of phase modulators
% Vapp - Drive voltage applied to modulator, V
%
% OUTPUT
% u - pulse at output

% uWave mod frequency, 2*pi*hz
omega_m = 2*pi * freq_m;    

% Vpi of phase modulator, V
V_pi_PM = 3;
V_pi_IM = 3;

% Convert power to voltage
V_IM = 2 * 10^((P_IM-10)/20);
V_PM = 2 * 10^((P_PM-10)/20);

% Phase set point of IM, radians
phi_dc = pi/2;  

% Drive Voltage of IM in units of Vpi
%vpi_ratio = 1/2;

% Number of intensity and phase modulators
N_IM = 1;
N_PM = 2;

% Intensity modulator pulse
E_IM = 1 + exp(1i * phi_dc + 1i * pi * V_IM / V_pi_IM .* cos(omega_m .* t));

% Phase modulator pulse
E_PM = exp(1i .* pi * V_PM / V_pi_PM .* cos(omega_m .* t));

% Combine the pulses
u = E_IM.^N_IM .* E_PM.^N_PM;

end % End function
