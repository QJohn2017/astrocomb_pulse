function [peak, fwhm] = prop_eom_comb(len_smf, edfa_power, len_hnlf, len_comp)

    %% Initialize Parameters

    freq_m = 16E9;                  % Modulator drive freq in Hz
    n_pulses = 2;                   % Number of pulses in the time window
    T = n_pulses / freq_m * 1E12;   % Time window, unit ps
    nt = 2^13;                      % Number of points
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
        'n_pulses', 't', 'f', 'lambda0', 'freq_m', 'nt', 'dt', 'n1', 'nplot', 'nz');    

    lambda_mask  = find(lambda > 0);  % find where the lambda > 0 for ploting

    %% Zeroth Step: EOM Comb

    % RF Power applied to Modulators, dBm
    P_PM = 25.5;
    P_IM = 7.5;

    % Number of intensity and phase modulators
    N_IM = 1;
    N_PM = 2;

    % Initialize EOM Pulse
    u_ini = EOM_pulse(t, freq_m*1E-12, P_IM, P_PM, N_IM, N_PM);

    %% First Step: SMF

    % Input average power, W
    avg_power = 0.001;

    % Length of fiber, m
    z = len_smf;

    [smf_t, smf_f, smf_power_t, smf_power_f, smf_phase_t, smf_phase_f, smf_fwhm] = ...
        prop_smf(u_ini, avg_power, z);

    %% Second Step: EOM Through EDFA. Note, accounts for 2 stage EDFA

    % Average power in each EDFA stage, W
    avg_power = [0.08, edfa_power];

    [edfa_t, edfa_f, edfa_power_t, edfa_power_f, edfa_phase_t, edfa_phase_f, edfa_fwhm] = ...
        prop_edfa(smf_t(:, end), avg_power);


    %% Third Step: EOM Through HNLF

    % Length of fiber, m
    z = len_hnlf;

    % Use the average power from the EDFA second stage
    avg_power = avg_power(2);

    [hnlf_t, hnlf_f, hnlf_power_t, hnlf_power_f, hnlf_phase_t, hnlf_phase_f, hnlf_fwhm] = ...
        prop_hnlf(edfa_t(:, end), avg_power, z);

    %% Fourth Step: Pulse Compressor

    % Length of fiber, m
    z = len_comp;
    
    % Introduce some fiber losses
    avg_power = 0.8 * avg_power;

    [comp_t, comp_f, comp_power_t, comp_power_f, comp_phase_t, comp_phase_f, comp_fwhm] = ...
        prop_smf(hnlf_t(:, end), avg_power, z);   
    
    peak = max(comp_power_t(:, end));
    fwhm = comp_fwhm;
    
end