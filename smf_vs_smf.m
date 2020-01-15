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

smf1s = [];
smf2s = [];
peaks = [];
fwhms = [];

for smf1 = 160:2.0:240
    for smf2 = 0.5:0.1:2.0
        [peak, fwhm] = prop_eom_comb(smf1, 3.5, 10, smf2);
        
        peaks = [peaks, peak];
        fwhms = [fwhms, fwhm];
        smf1s = [smf1s, smf1];
        smf2s = [smf2s, smf2];
    
        fprintf('SMF1 %.1f W, SMF2 %.1f m\n', smf1, smf2);
        fprintf('Compressor Pulse Width %.2f ps\n', fwhm);
        fprintf('Compressor Peak Power %.2f pJ\n', peak);
    end
end

save('smf_vs_smf_data.mat', 'peaks', 'fwhms', 'smf1s', 'smf2s');
