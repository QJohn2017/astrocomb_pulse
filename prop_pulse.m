function [out_t, out_f, power_t, power_f, phase_t, phase_f, fwhm] = ...
    prop_pulse(init_pulse_t, alpha, betap, gamma, avg_power, z)
    % t = Time domain
    % f = Frequency domain
    
    % Load global variables
    load('pulse_globalvars.mat',...
        'n_pulses', 't', 'f', 'lambda0', 'freq_m', 'nt', 'dt', 'n1', 'nplot', 'nz');	
    
    % Stepsize
    dz = z/nz;	
    % Vector of Points to plot (for heatmaps)
    zv = (z/nplot) * (0:nplot);	
    
    % Preallocate arrays
    out_t = zeros(length(t), length(zv));	
    out_f = zeros(length(t), length(zv));
    
    % Noise for Comb Frequency Spectrum
    A_f = randn(nt, 1) .* exp(-1i * randn(nt, 1));	
    
    % Noise for Comb Pulse
    A_t = fftshift(ifft(A_f)) * 1E-4;

    % Initial Time Pulse
    pulse_t = init_pulse_t + A_t;
    
    % Convert average power to energy, pJ
    energy = n_pulses * avg_power / freq_m * 1E12;
    
    % Initialize Time Pulse Array
    out_t(:, 1) = pulse_t .* sqrt(energy / trapz(t, abs(pulse_t).^2));
    
    % Initial Frequency Pulse
    pulse_f = fftshift(fft(out_t(:, 1)));
    
    % Initialize Frequency Pulse array, sqrt(pJ/THz)
    out_f(:, 1) = pulse_f .* sqrt(energy / trapz(f, abs(pulse_f).^2));
    
    % Propagate the pulse using ssprop
    for ii = 1:nplot  
      out_t(:, ii+1) = ssprop(out_t(:, ii), dt, dz, n1, alpha, betap, gamma);
      pulse_f = fftshift((fft(out_t(:, ii+1))));
      out_f(:, ii+1) = pulse_f .* sqrt(trapz(t, abs(out_t(:, ii+1)).^2) ./ trapz(f, abs(pulse_f).^2));
    end

    % Convert power to Friendly Units
    power_t = abs(out_t).^2;  % W
    power_sf = abs(out_f).^2 * freq_m;  % pW/THz
    power_slam = 3E8 ./ (lambda0.^2) .* power_sf*1E-15;  % W/nm
    power_f = 10 * log10(power_slam*1E3);  % dBm/nm
    
    % Calculate the phase
    phase_t = angle(out_t);
    phase_f = angle(out_f);
    
    % Fit a gaussian to the intensity pulse
    fit_mask = find((t > -5) & (t < 5));  % Only fit part of the waveform
    fit_t = fit(t(fit_mask), power_t(fit_mask, end), 'gauss1');
    fwhm = 2 * fit_t.c1;
end
