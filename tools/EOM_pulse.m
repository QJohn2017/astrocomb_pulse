function u = EOM_pulse(t, f_mod, Vapp)
%EOM_PULSE Summary of this function goes here
%   Equation taken from https://www.osapublishing.org/ol/abstract.cfm?uri=ol-35-19-3234
omega_m = 2*pi*f_mod;    %uWave mod frequency, 2*pi*hz
%Vapp = 15;   %RF voltage applied to Phase Modulator, V
V_pi = 3;   %Vpi of phase modulator, V
phi_dc = pi/2;  %Phase set point of AM, radians
vpi_ratio = 1/2; %1/2; %Drive Voltage of AM in units of Vpi
%vpi_ratio = 2/3.9;
N_IM = 1;
N_PM = 2;

E_IM = (1 + exp(1i*phi_dc).*exp(1i.*(pi*vpi_ratio).*cos(omega_m.*t))).^N_IM;    %Use 1 - exp... to center pulse
E_PM = (exp(1i.*(pi*Vapp/V_pi).*cos(omega_m.*t))).^N_PM;
u = E_IM.*E_PM;
end

