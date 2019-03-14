%% Pulse Propagation of an EOM Comb through Multiple Elements
% Originally Based off of Xianwen's code edited by Alex

% This program propagates an EOM comb through 4 parts: SMF, 2 stage EDFA (assumes
% anomolous dispersion in EDFA), HNLF, and pulse compressor

% For each element, you should set the dispersion (beta), loss (alpha),
% nonlinearity (gamma), and length. Note that the EOM comb and compresor
% have some other parameters as well.

clc
clear
close all
%load seaborn_colors.mat
current_directory = pwd;
addpath(fullfile(current_directory, 'tools'));   % path for Tools Directory
%addpath ('C:\Users\Alex\Documents\MATLAB\Astro_EOM_comb\Pulse Simulation\ssprop-3.0.1\tools')   % path for Tools Directory
% Don't need to add that path if this file, ssprop.m and EOM_pulse.m are on
% the same path already
%% Initialize Parameters

uWave_freq = 16E9;           %uWave drive freq in Hz
T = 125;                     % Time window, unit ps, should be long enough for ~2 cycles
nt = 2^14;                   % Number of points
dt = T/(nt-1);               % Timestep
t = ((1:nt)'-(nt+1)/2)*dt;   % Time vector  
w = wspace(T,nt);            % Angular frequency vector  
vs = fftshift(w/(2*pi));     % Shifted frequency (shifted for plotting), THz
lambda0 = 1550;              % Center Wavelength, nm
C = 3E8;                     % Speed of light, m/s
v0 = C/lambda0*1E-3;         % Center Frequency, THz
f = v0 + vs;                 % Shifted Center Frequency, THz
lambda=3E5./f;               % Shifted Center Wavelength,  nm
nz = 1000;                   % Number of calculation steps
nplot = 2;					 % Number of plots to make (only useful for heatmaps)
n1 = round(nz/nplot);		 % Number of steps per plot
save('pulse_globalvars.mat', ...
    't', 'f', 'lambda', 'uWave_freq', 'nt', 'dt', 'n1', 'nplot', 'nz');    %Save global Variables so they can be passed to the function


%% First Step: EOM Comb and SMF
z = 180;    % Total propogation distance, unit m
betap = 1.0*[0, 0, -0.0217, 5.0308e-05];     % Dispersion Beta Params as [0, 0, beta2, beta3] unit [0 0 ps2/m ps3/m] 
alpha = 0.18/1000/4.34;     % Loss coefficient, unit m-1. Leading number is loss in dB/km
ave_power = 0.001;      % Avg Power, W
P_RF = 25.5;      % RF Power applied to Phase Modulators
Vapp = (10^((P_RF-10)/20))*2;	%Drive voltage applied to modulator, V
Energy = ave_power/uWave_freq*1E12;	% Pulse Energy, pJ
fwhm = (1/uWave_freq)*1E12;	% Pulse width FWHM, ps. Defined as 50% of applied uWave rep rate
peak_power = (Energy*1E-12)/(fwhm*1E-12);	%Peak power in J(?)
gamma = 0.78*1E-3;	%Nonlinaer Coefficient of material Taken from https://ieeexplore.ieee.org/document/7764544

u_ini = sqrt(peak_power)*EOM_pulse(t, uWave_freq*1E-12, Vapp);  %Initialize EOM Pulse
[~, min_idx] = max(u_ini);  %Index of center of pulse
u_ini = circshift(u_ini, min_idx);  %Shift the waveform to the center of the time spectrum

tic;
disp('Propagating through SMF')
[smf_t, smf_f, smf_power_f, smf_power_t, smf_phi_t, smf_phi_f] = ...
    prop_pulse(u_ini, alpha, betap, gamma, Energy, z);
toc;

fit_mask = find((t > -5) & (t < 5));    %Only fit part of the waveform
smf_fit = fit(t(fit_mask), smf_power_t(fit_mask, end), 'gauss1');   %Gaussian fit of intensity

%% Second Step: EOM Through EDFA. Note, accounts for 2 stage EDFA

z = [8 10];	%Lenght of fiber in each step. Basically a guess
betap = [0, 0, 0, 5.0308e-05; 0, 0, -0.0217, 5.0308e-05];	% Dispersion Beta Params as [0, 0, beta2, beta3] unit [0 0 ps2/m ps3/m] 
alpha = 0.18/1000/4.34;	% Loss coefficient, unit m-1. Leading number is loss in dB/km
%alpha = (-35/z)/1000/4.34;
ave_power = [0.2, 2.5];	% Avg Power, W
fwhm = 2*smf_fit.c1;	% Pulse width FWHM, ps. Inherited from SMF
n2 = 2.5E-20;	% n2 of materialm m2/W
Aeff = pi*((6.5/2)^2)*1E-6*1E-6;	% Effective Area of mode, m^2
gamma = (n2/Aeff)*(2*pi/(lambda0*1E-9));    %nonlinear coefficient 1/(W m)
%u_ini = sqrt(peak_power)*smf_t(:, end);  %Initial Pulse from previous pulse

tic;
disp('Propagating through EDFA');
for i = 1:2     %loop twice becuase of 2 stage edfa
    Energy = ave_power(i)/uWave_freq*1E12; % Pulse Energy, pJ
    peak_power = (Energy*1E-12)/(fwhm*1E-12);   %Peak power in J(?)
    if i == 1
        u_ini = sqrt(peak_power)*smf_t(:, end);  %Initial Pulse from previous pulse
    else
        u_ini = sqrt(peak_power)*edfa_t(:, end);
    end
    [edfa_t, edfa_f, edfa_power_f, edfa_power_t, edfa_phi_t, edfa_phi_f] = ...
        prop_pulse(u_ini, alpha, betap(i, :), gamma, Energy, z(i));
end
toc;
edfa_fit = fit(t(fit_mask), edfa_power_t(fit_mask, end), 'gauss1'); %Fit of intensity pulse

%figure; plot(t, smf_power_t(:, end), t, edfa_power_t(:, end))

%% Third Step: EOM Through HNLF

z = 10;
alpha = 0.22/1000/4.34;       % Loss coefficient, unit m-1. Leading number is loss in dB/km
betap = [0, 0, 0.0023, -7.2317e-06];
gamma = 10.6*1E-3;
u_ini = edfa_t(:, end); %Initial pulse from previous pulse

tic;
disp('Propagating through HNLF');
[hnlf_t, hnlf_f, hnlf_power_f, hnlf_power_t, hnlf_phi_t, hnlf_phi_f] = ...
    prop_pulse(u_ini, alpha, betap, gamma, Energy, z);
toc;
hnlf_fit = fit(t(fit_mask), hnlf_power_t(fit_mask, end), 'gauss1'); %Fit of Intenisty FWHM

%figure; plot(t, edfa_power_t(:, end), t, hnlf_power_t(:, end));

%% Fourth Step: Pulse Compressor

z = 1.5;
betap = [0, 0, -0.0217, 5.0308e-05];     % Dispersion Beta Params as [0, 0, beta2, beta3] unit [0 0 ps2/m ps3/m] 
alpha = 0.18/1000/4.34;     % Loss coefficient, unit m-1. Leading number is loss in dB/km
gamma = 0.78*1E-3;	%Nonlinaer Coefficient of material Taken from https://ieeexplore.ieee.org/document/7764544
%Energy = ave_power(2)/uWave_freq*1E12; % Pulse Energy, pJ
%peak_power = (Energy*1E-12)/(2*hnlf_fit.c1*1E-12);   %Peak power in J(?)
u_ini = hnlf_t(:, end); %Initial pulse from previous pulse


tic;
disp('Propagating through Compressor');
[comp_t, comp_f, comp_power_f, comp_power_t, comp_phi_t, comp_phi_f] = ...
    prop_pulse(u_ini, alpha, betap, gamma, Energy, z);
toc;

comp_t = comp_t(:, end); comp_f = comp_f(:, end);
comp_power_f = comp_power_f(:, end); comp_power_t = comp_power_t(:, end);
comp_phi_t = comp_phi_t(:, end); comp_phi_f = comp_phi_f(:, end);
comp_fit = fit(t(fit_mask), comp_power_t(fit_mask, end), 'gauss1'); %Fit of Intenisty FWHM
fprintf('Final Pulse Width %.2f ps\n', 2*comp_fit.c1);

% beta_grating = -0.6; %Dispersion of Grating Compressor in ps^2 (note 1 ps^2 = 10^6 fs^2). Should be negative
% phi = exp(-1i*beta_grating.*(f-v0).^2);	%Phase added by grating compressor
% E_ini = 0.8*hnlf_t(:,end);  %Estimate of loss of compressor
% E_in = fftshift(fft(E_ini));	%Input is FFT of hnlf time spectrum
% E_comp = E_in.*phi;	%compensated frequency spectrum
% comp_t = ifft(E_comp);	%Time spectrum of compressor
% comp_power_t = abs(comp_t).^2;	%Temporal power spectrum of compressor
% comp_f = sqrt(trapz(t,abs(comp_t).^2)./(trapz(f,abs(E_comp).^2))).*E_comp;  %Frequency Spectrum of Compressor
% 
% comp_power_sf = abs(comp_f).^2*uWave_freq;   %Frequency Spectrum in pW/THz
% comp_power_slam = 3E8./(lambda.^2).*comp_power_sf*1E-15;  %Frequency Spectrum in W/nm
% comp_power_f = 10*log10(comp_power_slam*1E3); %Frequency Spectrum in dBm/nm
% comp_phi_t = angle(comp_t); %Phase of time spectrum
% comp_phi_f = angle(comp_f); %Phasse of frequency spectrum
% 
% comp_fit = fit(t(fit_mask), comp_power_t(fit_mask, end), 'gauss1'); 
% fprintf('Final Pulse Width %.2f ps\n', 2*comp_fit.c1);
%figure; plot(t, hnlf_power_t(:, end), t, comp_power_t); xlim([-1 1])

%% Plot Everything
xlim_lambda = [1400 1700];

figure
scrsz = get(groot,'ScreenSize');
set(gcf, 'position',scrsz);

subplot(2, 1, 1)    %Time Plot
    ax = plot(t, (abs(smf_power_t(:, 1)).^2)./max(abs(smf_power_t(:, 1)).^2), ...
        t, smf_power_t(:, end)./max(smf_power_t(:, end)), ...
        t, edfa_power_t(:, end)./max(edfa_power_t(:, end)), ...
        t, hnlf_power_t(:, end)./max(hnlf_power_t(:, end)), ...
        t, comp_power_t./max(comp_power_t));
    xlabel('Time (ps)')
    ylabel('Norm Power')
    xlim([-3 3]); ylim([0 1.5]);
    title('Normalized Pulse Intensity')
    legend('Initial EOM Pulse', 'After SMF', 'After EDFA', 'After HNLF', 'After Grating Compressor');
    set(ax, 'LineWidth', 2)
    
subplot(2, 1, 2)    %Frequency Spectrum
    plot(lambda, smf_power_f(:, 1), ...
        lambda, smf_power_f(:, end), ...
        lambda, edfa_power_f(:, end), ...
        lambda, hnlf_power_f(:, end), ...
        lambda, comp_power_f);
    xlabel('Wavelength (nm)');
    ylabel('Power Density (dBm/nm)');
    xlim(xlim_lambda)
    title('Comb Spectrum')
    
%% Time Intensity and Phase
[~, ref_idx] = min(abs(t-0));   %Reference for phase angle

figure
scrsz = get(groot,'ScreenSize');
set(gcf, 'position',scrsz);

subplot(5, 1, 1)
    plot_yy(t, smf_power_t(:, 1), unwrap(smf_phi_t(:, 1) - smf_phi_t(ref_idx, 1)))
    title(sprintf('Initial EOM Comb, FWHM = %.1f ps', (1/(uWave_freq*2))*1E12)); xlim([-5 5])
subplot(5, 1, 2)
    plot_yy(t, smf_power_t(:, end), unwrap(smf_phi_t(:, end) - smf_phi_t(ref_idx, end)))
    title(sprintf('After SMF, FWHM = %.1f ps', 2*smf_fit.c1)); xlim([-5 5])
subplot(5, 1 , 3)
    plot_yy(t, edfa_power_t(:, end), unwrap(edfa_phi_t(:, end) - edfa_phi_t(ref_idx, end)))
    title(sprintf('After EDFA, FWHM = %.1f ps', 2*edfa_fit.c1)); xlim([-5 5])
subplot(5, 1, 4)
    plot_yy(t, hnlf_power_t(:, end), unwrap(hnlf_phi_t(:, end) - hnlf_phi_t(ref_idx, end)));
    title(sprintf('After HNLF, FWHM = %.1f ps', 2*hnlf_fit.c1)); xlim([-5 5])
subplot(5, 1, 5)
    plot_yy(t, comp_power_t(:, end), unwrap(comp_phi_t(:, end) - comp_phi_t(ref_idx, end)));
    title(sprintf('After Compressor, FWHM = %.2f ps', 2*comp_fit.c1)); xlim([-5 5])
    
%% Spectral Intensity and Phase
%[~, ref_idx] = min(abs(lambda-lambda0));

figure
scrsz = get(groot,'ScreenSize');
set(gcf, 'position',scrsz);

subplot(5, 1, 1)
    plot_yy(lambda, smf_power_f(:, 1), unwrap(smf_phi_f(:, 1) - smf_phi_f(ref_idx, 1)))
    title('Initial EOM Comb'); xlim(xlim_lambda)
subplot(5, 1, 2)
    plot_yy(lambda, smf_power_f(:, end), unwrap(smf_phi_f(:, end)- smf_phi_f(ref_idx, end)))
    title('After SMF'); xlim(xlim_lambda)
subplot(5, 1 , 3)
    plot_yy(lambda, edfa_power_f(:, end), unwrap(edfa_phi_f(:, end) - edfa_phi_f(ref_idx, end)))
    title('After EDFA'); xlim(xlim_lambda)
subplot(5, 1, 4)
    plot_yy(lambda, hnlf_power_f(:, end), unwrap(hnlf_phi_f(:, end) - hnlf_phi_f(ref_idx, end)));
    title('After HNLF'); xlim(xlim_lambda)
subplot(5, 1, 5)
    plot_yy(lambda, comp_power_f(:, end), unwrap(comp_phi_f(:, end) - comp_phi_f(ref_idx, end)));
    title('After Compressor'); xlim(xlim_lambda)
%% Functions

function [u_out, U_out, power_slam_dB, power_t, temporal_phase, spectral_phase] = prop_pulse(...
    u_ini, alpha, betap, gamma, Energy, z)
    % u = Time domain
    % U = Freq domain
    load('pulse_globalvars.mat');	%Load global variables
    dz = z/nz;	% Stepsize
    zv = (z/nplot)*(0:nplot);	% Vector of Points to plot (for heatmaps)
    u_out = zeros(length(t),length(zv));	% Preallocate arrays
    U_out = zeros(length(t),length(zv));
    

    A_f = randn(nt,1).*exp(-1i*randn(nt,1));	%Noise for Comb Frequency Spectrum
    A_t = fftshift(ifft(A_f))*1E-4;	%Initial Comb Spectrum

    u0 = u_ini + A_t; %Initial Time Pulse
    u_out(:,1) = sqrt(Energy/trapz(t,abs(u0).^2)).*u0;	%Initialize Time Pulse Array
    U0 = fftshift(fft(u_out(:,1)));	%Initial Frequency Pulse
    U_out(:,1) = sqrt(Energy/(trapz(f,abs(U0).^2))).*U0;	%Initialize Frequency Pulse array, sqrt(pJ/THz)
    for ii = 1:nplot  
      u_out(:,ii+1) = ssprop(u_out(:,ii),dt,dz,n1,alpha,betap,gamma);
        % %     %U(:,ii+1) = fftshift(abs(dt*fft(u(:,ii+1))/sqrt(2*pi)).^2);
        % % %    U(:,ii+1) = fftshift(abs(fft(u(:,ii+1))/nt).^2);
      U0 = fftshift((fft(u_out(:,ii+1))));
      U_out(:,ii+1) = sqrt(trapz(t,abs(u_out(:,ii+1)).^2)./(trapz(f,abs(U0).^2))).*U0;
    end

    % Convert to Friendly Units
     power_sf = abs(U_out).^2*uWave_freq;   % pW/THz
     power_slam = 3E8./(lambda.^2).*power_sf*1E-15;  %W/nm
     power_slam_dB = 10*log10(power_slam*1E3); %dBm/nm
     power_t=abs(u_out).^2;   %W
     temporal_phase = angle(u_out);
     spectral_phase = angle(U_out);
end

function plot_yy(xdata, ydata1, ydata2)
    yyaxis left
    plot(xdata, ydata1, '-');
    ylabel('Intensity (au)')
    yyaxis right
    plot(xdata, ydata2, '--');
    ylabel('Phase (rad)');
end