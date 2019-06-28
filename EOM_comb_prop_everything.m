%% Pulse Propagation of an EOM Comb through Multiple Elements
% Originally Based off of Xianwen's code edited by Alex and Ryan

% This program propagates an EOM comb through 4 parts: SMF, 2 stage EDFA (assumes
% anomolous dispersion in EDFA), HNLF, and pulse compressor

% For each element, you should set the dispersion (beta), loss (alpha),
% nonlinearity (gamma), and length. Note that the EOM comb and compresor
% have some other parameters as well.

clc
clear
close all
cur_dir = pwd;
% path for Tools Directory (ssprop.m and EOM_pulse.m)
addpath(fullfile(cur_dir, 'tools'));

%% Initialize Parameters

freq_m = 16E9;                  % Modulator drive freq in Hz
T = 2.0 / freq_m * 1E12;        % Time window, unit ps, should be long enough for ~2 cycles
nt = 2^16;                      % Number of points
dt = T / (nt-1);                % Timestep
t = ((1:nt)' - (nt+1)/2) * dt;  % Time vector  
w = wspace(T, nt);              % Angular frequency vector  
vs = fftshift(w / (2*pi));      % Shifted frequency (shifted for plotting), THz
lambda0 = 1550;                 % Center Wavelength, nm
C = 3E8;                        % Speed of light, m/s
v0 = C / lambda0*1E-3;          % Center Frequency, THz
f = v0 + vs;                    % Shifted Center Frequency, THz
lambda = 3E5 ./ f;              % Shifted Center Wavelength,  nm
nz = 1000;                      % Number of calculation steps
nplot = 2;					    % Number of plots to make (only useful for heatmaps)
n1 = round(nz/nplot);		    % Number of steps per plot

% Save global Variables so they can be passed to the function
save('pulse_globalvars.mat', ...
    't', 'f', 'lambda', 'freq_m', 'nt', 'dt', 'n1', 'nplot', 'nz');    

lambda_mask  = find(lambda > 0);  % find where the lambda > 0 for ploting

%% Zeroth Step: EOM Comb

% Avg Power, W
avg_power = 0.001;
% RF Power applied to Modulators
P_PM = 25.5;
P_IM = 7.5;
% Pulse energy, pJ
energy = avg_power / freq_m * 1E12;
% Pulse width FWHM, ps. Defined as 50% of applied uWave rep rate
fwhm = (1 / freq_m) * 1E12;
% Peak power in W
peak_power = (energy*1E-12) / (fwhm*1E-12);

% Initialize EOM Pulse
u_ini = sqrt(peak_power) * EOM_pulse(t, freq_m*1E-12, P_IM, P_PM);

% Index of center of pulse
[~, min_idx] = max(u_ini);

% Shift the waveform to the center of the time spectrum
u_ini = circshift(u_ini, min_idx);  
%u_ini = circshift(u_ini, length(t)/2);


%% First Step: SMF

% Length of fiber, unit m
z = 180;

[smf_t, smf_f, smf_power_t, smf_power_f, smf_phase_t, smf_phase_f, smf_fwhm] = ...
    prop_smf(u_ini, avg_power, z);

fprintf('SMF Pulse Width %.2f ps\n', smf_fwhm);

%% Second Step: EOM Through EDFA. Note, accounts for 2 stage EDFA

% Average power after EDFA second stage, W
avg_power = 3.5;

[edfa_t, edfa_f, edfa_power_t, edfa_power_f, edfa_phase_t, edfa_phase_f, edfa_fwhm] = ...
    prop_edfa(smf_t(:, end), avg_power, lambda0);

fprintf('EDFA Pulse Width %.2f ps\n', edfa_fwhm);

%% Third Step: EOM Through HNLF

% Length of fiber, unit m
z = 10;

[hnlf_t, hnlf_f, hnlf_power_t, hnlf_power_f, hnlf_phase_t, hnlf_phase_f, hnlf_fwhm] = ...
    prop_hnlf(edfa_t(:, end), avg_power, z);

%% Fourth Step: Pulse Compressor

% Length of fiber, unit m
z = 1.0;

[comp_t, comp_f, comp_power_t, comp_power_f, comp_phase_t, comp_phase_f, comp_fwhm] = ...
    prop_smf(hnlf_t(:, end), avg_power, z);

fprintf('Compressor Pulse Width %.2f ps\n', comp_fwhm);
fprintf('Compressor Peak Power %.2f pJ\n', max(comp_power_t(:, end)));

%% Fifth Step: Waveguide

% Length of waveguide, unit m
z = 10e-3;

[wg_t, wg_f, wg_power_t, wg_power_f, wg_phase_t, wg_phase_f, wg_fwhm] = ...
    prop_waveguide(comp_t(:, end), avg_power, z, lambda0);

fprintf('Final Pulse Width %.2f ps\n', wg_fwhm);

% % Pulse energy, pJ
% energy = fiber_loss * sqrt(trans) * avg_power / freq_m*1E12;
% 
% % Peak power in J(?)
% peak_power = (energy*1E-12) / (comp_fwhm*1E-12); 
% fprintf('Peak Power into Waveguide %.2f J\n', peak_power);

%% Plot Everything
xlim_lambda = [1400 1700];

figure
scrsz = get(groot,'ScreenSize');
set(gcf, 'position',scrsz);

% Time Plot
subplot(2, 1, 1)    
    ax = plot(t, (abs(smf_power_t(:, 1)).^2) ./ max(abs(smf_power_t(:, 1)).^2), ...
        t, smf_power_t(:, end) ./ max(smf_power_t(:, end)), ...
        t, edfa_power_t(:, end) ./ max(edfa_power_t(:, end)), ...
        t, hnlf_power_t(:, end) ./ max(hnlf_power_t(:, end)), ...
        t, comp_power_t(:, end) ./ max(comp_power_t(:, end)), ...
        t, wg_power_t(:, end) ./ max(wg_power_t(:, end)));
    xlabel('Time (ps)')
    ylabel('Norm Power')
    xlim([-3 3]); ylim([0 1.5]);
    title('Normalized Pulse Intensity')
    legend('Initial EOM Pulse', 'After SMF', 'After EDFA', 'After HNLF', 'After Grating Compressor');
    set(ax, 'LineWidth', 2)
    
% Frequency Spectrum
subplot(2, 1, 2)   
    plot(lambda, smf_power_f(:, 1), ...
        lambda, smf_power_f(:, end), ...
        lambda, edfa_power_f(:, end), ...
        lambda, hnlf_power_f(:, end), ...
        lambda, comp_power_f(:, end), ...
        lambda, wg_power_f(:, end));
    xlabel('Wavelength (nm)');
    ylabel('Power Density (dBm/nm)');
    xlim(xlim_lambda)
    title('Comb Spectrum')
    
%% Time Intensity and Phase
[~, ref_idx] = min(abs(t-0));   %Reference for phase angle

figure
scrsz = get(groot,'ScreenSize');
set(gcf, 'position',scrsz);

subplot(6, 1, 1)
    plot_yy(t, smf_power_t(:, 1), unwrap(smf_phase_t(:, 1) - smf_phase_t(ref_idx, 1)))
    title(sprintf('Initial EOM Comb, FWHM = %.1f ps', (1/(freq_m*2))*1E12)); xlim([-10 10])
subplot(6, 1, 2)
    plot_yy(t, smf_power_t(:, end), unwrap(smf_phase_t(:, end) - smf_phase_t(ref_idx, end)))
    title(sprintf('After SMF, FWHM = %.1f ps', smf_fwhm)); xlim([-5 5])
subplot(6, 1 , 3)
    plot_yy(t, edfa_power_t(:, end), unwrap(edfa_phase_t(:, end) - edfa_phase_t(ref_idx, end)))
    title(sprintf('After EDFA, FWHM = %.1f ps', edfa_fwhm)); xlim([-5 5])
subplot(6, 1, 4)
    plot_yy(t, hnlf_power_t(:, end), unwrap(hnlf_phase_t(:, end) - hnlf_phase_t(ref_idx, end)));
    title(sprintf('After HNLF, FWHM = %.1f ps', hnlf_fwhm)); xlim([-5 5])
subplot(6, 1, 5)
    plot_yy(t, comp_power_t(:, end), unwrap(comp_phase_t(:, end) - comp_phase_t(ref_idx, end)));
    title(sprintf('After Compressor, FWHM = %.2f ps', comp_fwhm)); xlim([-5 5])
subplot(6, 1, 6)
    plot_yy(t, wg_power_t(:, end), unwrap(wg_phase_t(:, end) - wg_phase_t(ref_idx, end)));
    title(sprintf('After Waveguide, FWHM = %.2f ps', wg_fwhm)); xlim([-5 5])
    
%% Spectral Intensity and Phase
%[~, ref_idx] = min(abs(lambda-lambda0));

figure
scrsz = get(groot,'ScreenSize');
set(gcf, 'position',scrsz);

subplot(6, 1, 1)
    plot_yy(lambda, smf_power_f(:, 1), unwrap(smf_phase_f(:, 1) - smf_phase_f(ref_idx, 1)))
    title('Initial EOM Comb'); xlim(xlim_lambda)
subplot(6, 1, 2)
    plot_yy(lambda, smf_power_f(:, end), unwrap(smf_phase_f(:, end)- smf_phase_f(ref_idx, end)))
    title('After SMF'); xlim(xlim_lambda)
subplot(6, 1 , 3)
    plot_yy(lambda, edfa_power_f(:, end), unwrap(edfa_phase_f(:, end) - edfa_phase_f(ref_idx, end)))
    title('After EDFA'); xlim(xlim_lambda)
subplot(6, 1, 4)
    plot_yy(lambda, hnlf_power_f(:, end), unwrap(hnlf_phase_f(:, end) - hnlf_phase_f(ref_idx, end)));
    title('After HNLF'); xlim(xlim_lambda)
subplot(6, 1, 5)
    plot_yy(lambda, comp_power_f(:, end), unwrap(comp_phase_f(:, end) - comp_phase_f(ref_idx, end)));
    title('After Compressor'); xlim(xlim_lambda)    
subplot(6, 1, 6)
    plot_yy(lambda, wg_power_f(:, end), unwrap(wg_phase_f(:, end) - wg_phase_f(ref_idx, end)));
    title('After Waveguide'); xlim([500 2000])
    
%% Functions

function plot_yy(xdata, ydata1, ydata2)
    yyaxis left
    plot(xdata, ydata1, '-');
    ylabel('Intensity (au)')
    yyaxis right
    plot(xdata, ydata2, '--');
    ylabel('Phase (rad)');
end
