function pulse = EOM_pulse(t, freq_m, P_IM, P_PM, N_IM, N_PM)
% Produce a pulse from electro-optic intensity and phase modulators
% Equation taken from https://www.osapublishing.org/ol/abstract.cfm?uri=ol-35-19-3234
%
% INPUT
%
% t - time series, ps
% freq_m - drive frequency, THz
% P_IM - drive power of intensity modulators, dBm
% P_PM - drive power of phase modulators, dBm
% N_IM - number of intensity modulators
% N_PM - number of phase modulators
%
% OUTPUT
% 
% u - pulse at output

% Vpi of modulators, V
V_pi_PM = 3;
V_pi_IM = 3;

% Convert power to voltage
V_IM = 2 * 10^((P_IM-10)/20);
V_PM = 2 * 10^((P_PM-10)/20);

% Phase set point of IM, radians
phi_dc = pi/2;  

% Drive Voltage of IM in units of Vpi
%vpi_ratio = 1/2;

% Intensity modulator pulse
E_IM = 1 + exp(1i*phi_dc + 1i*pi * V_IM/V_pi_IM .* cos(2*pi*freq_m .* t));

% Phase modulator pulse
E_PM = exp(1i*pi * V_PM/V_pi_PM .* cos(2*pi*freq_m .* t));

% Combine the pulses
pulse = E_IM.^N_IM .* E_PM.^N_PM;

% Index of center of pulse
[~, max_idx] = max(abs(pulse));

% Shift the waveform to the center of the time spectrum
pulse = circshift(pulse, length(t)/2 - max_idx);

end % End function
