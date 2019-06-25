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
ave_power = 0.001;
% RF Power applied to Modulators
P_PM = 25.5;
P_IM = 7.5;
% Pulse energy, pJ
energy = ave_power / freq_m * 1E12;
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

% Total propogation distance, unit m
z = 200;
% Dispersion Beta Params as [0, 0, beta2, beta3] unit [0 0 ps2/m ps3/m] 
betap = 1.0 * [0, 0, -0.0217, 5.0308e-05];
% Loss coefficient, unit m-1. Leading number is loss in dB/km
alpha = 0.18 / 1000 / 4.34;
% Nonlinear coefficient of material taken from https://ieeexplore.ieee.org/document/7764544
gamma = 0.78*1E-3;

disp('Propagating through SMF')
tic;
[smf_t, smf_f, smf_power_f, smf_power_t, smf_phi_t, smf_phi_f] = ...
    prop_pulse(u_ini, alpha, betap, gamma, energy, z);
toc;

% Gaussian fit of intensity
fit_mask = find((t > -5) & (t < 5));  % Only fit part of the waveform
smf_fit = fit(t(fit_mask), smf_power_t(fit_mask, end), 'gauss1');
fprintf('SMF Pulse Width %.2f ps\n', 2*smf_fit.c1); 

%% Second Step: EOM Through EDFA. Note, accounts for 2 stage EDFA

% Length of fiber (m) in each step. Basically a guess
z = [8 10];
% Dispersion Beta Params as [0, 0, beta2, beta3] unit [0 0 ps2/m ps3/m] 
betap = [0, 0, 0, 5.0308e-05;
         0, 0, -0.0217, 5.0308e-05];
% Loss coefficient, unit m-1. Leading number is loss in dB/km
alpha = 0.18 / 1000 / 4.34;
%alpha = (-35/z)/1000/4.34;
% Avg Power, W
ave_power = [0.2, 3.75];
% Pulse width FWHM, ps. Inherited from SMF
fwhm = 2*smf_fit.c1;
% n2 of material m2/W
n2 = 2.5E-20;
% Effective Area of mode, m^2
Aeff = pi*((6.5/2)^2)*1E-6*1E-6;
% Nonlinear coefficient 1/(W m)
gamma = (n2/Aeff) * (2*pi / (lambda0*1E-9));

disp('Propagating through EDFA');
tic;
% Loop twice becuase of 2 stage edfa
for i = 1:2     
    % Pulse energy, pJ
    energy = ave_power(i) / freq_m*1E12; 
    % Peak power in J(?)
    peak_power = (energy*1E-12) / (fwhm*1E-12);
    
    % Initialize pulse from previous pulse
    if i == 1
        u_ini = sqrt(peak_power) * smf_t(:, end);  
    else
        u_ini = sqrt(peak_power) * edfa_t(:, end);
    end
    
    [edfa_t, edfa_f, edfa_power_f, edfa_power_t, edfa_phi_t, edfa_phi_f] = ...
        prop_pulse(u_ini, alpha, betap(i, :), gamma, energy, z(i));
end
toc;

% Fit of intensity pulse
edfa_fit = fit(t(fit_mask), edfa_power_t(fit_mask, end), 'gauss1');
fprintf('EDFA Pulse Width %.2f ps\n', 2*edfa_fit.c1); 

%% Third Step: EOM Through HNLF

% Length of fiber, unit m
z = 10;
% Loss coefficient, unit m-1. Leading number is loss in dB/km
alpha = 0.22 / 1000 / 4.34;
% Dispersion Beta Params as [0, 0, beta2, beta3] unit [0 0 ps2/m ps3/m] 
betap = [0, 0, 0.0023, -7.2317e-06];
% Nonlinear coefficient 1/(W m)
gamma = 10.6*1E-3;

% Initial pulse from previous pulse
u_ini = edfa_t(:, end);

disp('Propagating through HNLF');
tic;
[hnlf_t, hnlf_f, hnlf_power_f, hnlf_power_t, hnlf_phi_t, hnlf_phi_f] = ...
    prop_pulse(u_ini, alpha, betap, gamma, energy, z);
toc;

% Fit of Intenisty FWHM
hnlf_fit = fit(t(fit_mask), hnlf_power_t(fit_mask, end), 'gauss1');

%% Fourth Step: Pulse Compressor

% Length of fiber, unit m
z = 1.1;
% Dispersion Beta Params as [0, 0, beta2, beta3] unit [0 0 ps2/m ps3/m]
betap = [0, 0, -0.0217, 5.0308e-05];
% Loss coefficient, unit m-1. Leading number is loss in dB/km
alpha = 0.18/1000/4.34;
% Nonlinaer coefficient of material taken from https://ieeexplore.ieee.org/document/7764544
gamma = 0.78*1E-3;	
%energy = ave_power(2)/freq_m*1E12; % Pulse energy, pJ
%peak_power = (energy*1E-12)/(2*hnlf_fit.c1*1E-12);   %Peak power in J(?)

% Initialize pulse from previous pulse
u_ini = hnlf_t(:, end);

disp('Propagating through Compressor');
tic;
[comp_t, comp_f, comp_power_f, comp_power_t, comp_phi_t, comp_phi_f] = ...
    prop_pulse(u_ini, alpha, betap, gamma, energy, z);
toc;

comp_fit = fit(t(fit_mask), comp_power_t(fit_mask, end), 'gauss1'); %Fit of Intenisty FWHM
fprintf('Compressor Pulse Width %.2f ps\n', 2*comp_fit.c1);

%% Fifth Step: HNLF #2

% Length of fiber, unit m
z = 10;
% Loss coefficient, unit m-1. Leading number is loss in dB/km
alpha = 0.22 / 1000 / 4.34;
% Dispersion Beta Params as [0, 0, beta2, beta3] unit [0 0 ps2/m ps3/m] 
betap = [0, 0, 0.0023, -7.2317e-06];
% Nonlinear coefficient 1/(W m)
gamma = 10.6*1E-3;

% Initial pulse from previous pulse
u_ini = comp_t(:, end);

disp('Propagating through HNLF');
tic;
[hnlf2_t, hnlf2_f, hnlf2_power_f, hnlf2_power_t, hnlf2_phi_t, hnlf2_phi_f] = ...
    prop_pulse(u_ini, alpha, betap, gamma, energy, z);
toc;

% Fit of Intenisty FWHM
hnlf2_fit = fit(t(fit_mask), hnlf2_power_t(fit_mask, end), 'gauss1');

%% Plot Everything
xlim_lambda = [1400 1700];

figure
scrsz = get(groot,'ScreenSize');
set(gcf, 'position',scrsz);

% Time Plot
subplot(2, 1, 1)    
    ax = plot(t, (abs(smf_power_t(:, 1)).^2)./max(abs(smf_power_t(:, 1)).^2), ...
        t, smf_power_t(:, end)./max(smf_power_t(:, end)), ...
        t, edfa_power_t(:, end)./max(edfa_power_t(:, end)), ...
        t, hnlf_power_t(:, end)./max(hnlf_power_t(:, end)), ...
        t, comp_power_t(:, end)./max(comp_power_t(:, end)), ...
        t, hnlf2_power_t./max(hnlf2_power_t(:, end)));
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
        lambda, hnlf2_power_f(:, end));
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
    plot_yy(t, smf_power_t(:, 1), unwrap(smf_phi_t(:, 1) - smf_phi_t(ref_idx, 1)))
    title(sprintf('Initial EOM Comb, FWHM = %.1f ps', (1/(freq_m*2))*1E12)); xlim([-10 10])
subplot(6, 1, 2)
    plot_yy(t, smf_power_t(:, end), unwrap(smf_phi_t(:, end) - smf_phi_t(ref_idx, end)))
    title(sprintf('After SMF, FWHM = %.1f ps', 2*smf_fit.c1)); xlim([-5 5])
subplot(6, 1 , 3)
    plot_yy(t, edfa_power_t(:, end), unwrap(edfa_phi_t(:, end) - edfa_phi_t(ref_idx, end)))
    title(sprintf('After EDFA, FWHM = %.1f ps', 2*edfa_fit.c1)); xlim([-5 5])
subplot(6, 1, 4)
    plot_yy(t, hnlf_power_t(:, end), unwrap(hnlf_phi_t(:, end) - hnlf_phi_t(ref_idx, end)));
    title(sprintf('After HNLF, FWHM = %.1f ps', 2*hnlf_fit.c1)); xlim([-5 5])
subplot(6, 1, 5)
    plot_yy(t, comp_power_t(:, end), unwrap(comp_phi_t(:, end) - comp_phi_t(ref_idx, end)));
    title(sprintf('After Compressor, FWHM = %.2f ps', 2*comp_fit.c1)); xlim([-5 5])
subplot(6, 1, 6)
    plot_yy(t, hnlf2_power_t(:, end), unwrap(hnlf2_phi_t(:, end) - hnlf2_phi_t(ref_idx, end)));
    title(sprintf('After Waveguide, FWHM = %.2f ps', 2*hnlf2_fit.c1)); xlim([-5 5])
    
%% Spectral Intensity and Phase
%[~, ref_idx] = min(abs(lambda-lambda0));

figure
scrsz = get(groot,'ScreenSize');
set(gcf, 'position',scrsz);

subplot(6, 1, 1)
    plot_yy(lambda, smf_power_f(:, 1), unwrap(smf_phi_f(:, 1) - smf_phi_f(ref_idx, 1)))
    title('Initial EOM Comb'); xlim(xlim_lambda)
subplot(6, 1, 2)
    plot_yy(lambda, smf_power_f(:, end), unwrap(smf_phi_f(:, end)- smf_phi_f(ref_idx, end)))
    title('After SMF'); xlim(xlim_lambda)
subplot(6, 1 , 3)
    plot_yy(lambda, edfa_power_f(:, end), unwrap(edfa_phi_f(:, end) - edfa_phi_f(ref_idx, end)))
    title('After EDFA'); xlim(xlim_lambda)
subplot(6, 1, 4)
    plot_yy(lambda, hnlf_power_f(:, end), unwrap(hnlf_phi_f(:, end) - hnlf_phi_f(ref_idx, end)));
    title('After HNLF'); xlim(xlim_lambda)
subplot(6, 1, 5)
    plot_yy(lambda, comp_power_f(:, end), unwrap(comp_phi_f(:, end) - comp_phi_f(ref_idx, end)));
    title('After Compressor'); xlim(xlim_lambda)    
subplot(6, 1, 6)
    %plot_yy(lambda, hnlf2_power_f, unwrap(hnlf2_phi_f(:, end) - hnlf2_phi_f(ref_idx, end)));
    plot(lambda, hnlf2_power_f(:, end));
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
