function [out_t, out_f, power_f, power_t, temporal_phase, spectral_phase] = prop_pulse(...
    init_pulse_t, alpha, betap, gamma, energy, z)
    % t = Time domain
    % f = Freq domain
    
    % Load global variables
    load('pulse_globalvars.mat',...
        't', 'f', 'lambda', 'freq_m', 'nt', 'dt', 'n1', 'nplot', 'nz');	
    
    % Stepsize
    dz = z/nz;	
    % Vector of Points to plot (for heatmaps)
    zv = (z/nplot) * (0:nplot);	
    
    % Preallocate arrays
    out_t = zeros(length(t), length(zv));	
    out_f = zeros(length(t), length(zv));
    
    % Noise for Comb Frequency Spectrum
    A_f = randn(nt, 1) .* exp(-1i * randn(nt, 1));	
    
    % Initial Comb Spectrum
    A_t = fftshift(ifft(A_f))*1E-4;	

    % Initial Time Pulse
    pulse_t = init_pulse_t + A_t; 
    
    % Initialize Time Pulse Array
    out_t(:,1) = pulse_t .* sqrt(energy / trapz(t, abs(pulse_t).^2));
    
    % Initial Frequency Pulse
    pulse_f = fftshift(fft(out_t(:, 1)));
    % Initialize Frequency Pulse array, sqrt(pJ/THz)
    out_f(:,1) = pulse_f .* sqrt(energy / (trapz(f, abs(pulse_f).^2)));
    
    % Propagate the pulse using ssprop
    for ii = 1:nplot  
      out_t(:, ii+1) = ssprop(out_t(:, ii), dt, dz, n1, alpha, betap, gamma);
        % U(:,ii+1) = fftshift(abs(dt*fft(u(:,ii+1))/sqrt(2*pi)).^2);
        % U(:,ii+1) = fftshift(abs(fft(u(:,ii+1))/nt).^2);
      pulse_f = fftshift((fft(out_t(:, ii+1))));
      out_f(:, ii+1) = pulse_f .* sqrt(trapz(t, abs(out_t(:, ii+1)).^2) ./ trapz(f, abs(pulse_f).^2));
    end

    % Convert to Friendly Units
    power_sf = abs(out_f).^2 * freq_m;  % pW/THz
    power_slam = 3E8 ./ (lambda.^2) .* power_sf*1E-15;  % W/nm
    power_f = 10 * log10(power_slam*1E3);  % dBm/nm
    power_t = abs(out_t).^2;  % W
    temporal_phase = angle(out_t);
    spectral_phase = angle(out_f);
end
