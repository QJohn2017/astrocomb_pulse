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

%% Iterate through the EOM Comb with different parameters

powers = [];
smfs = [];
peaks = [];
fwhms = [];

for edfa_power = 2.5:0.1:5.0
    for len_smf = 0.5:0.1:2.0
        [peak, fwhm] = prop_eom_comb(200, edfa_power, 10, len_smf);
        
        peaks = [peaks, peak];
        fwhms = [fwhms, fwhm];
        powers = [powers, edfa_power];
        smfs = [smfs, len_smf];
    
        fprintf('EDFA %.1f W, SMF %.1f m\n', edfa_power, len_smf);
        fprintf('Compressor Pulse Width %.2f ps\n', fwhm);
        fprintf('Compressor Peak Power %.2f pJ\n', peak);
    end
end

save('smf_vs_edfa_data.mat', 'peaks', 'fwhms', 'powers', 'smfs');
