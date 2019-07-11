function [beta2, beta3] = D_to_beta(D, d_D)

% Convert dispersion (D) and dispersion slope (d_D) into the
% appropriate beta parameters (GVD and TOD)
% 
% Basic derivation:
% https://www.newport.com/n/the-effect-of-dispersion-on-ultrashort-pulses
%
% INPUT
%
% D - dispersion, ps/(nm km)
% d_D - dispersion slope, ps/(nm^2 km)
%
% OUTPUT
%
% beta2 - group velocity dispersion (GVD), ps/m
% beta3 - third order dispersion (TOD), ps^2/m

% Standard parameters
c = 3E8 * 1E9 / 1E12; % Speed of light, m/s -> nm/ps
lambda = 1550; % Center wavelength, nm

% Dispersion in ps/(nm km) -> ps/(nm m)
D = D / 1E3;

% Second beta parameter in ps^2/m
beta2 = - D * lambda^2 / (2*pi*c);

% Dispersion slope in ps/(nm^2 km) -> ps/(nm^2 m)
d_D = d_D / 1E3;

% Third beta parameter in ps^3/m
beta3 = lambda^4 / (2*pi*c)^2 * (d_D + 2*D/lambda);



% beta3 = lambda^4 / (2*pi*c)^2 * (3*D/lambda + lambda/(c*D) * d_D);

end

function [beta2, beta3] = D_to_beta2(D, d_D) %Input D (dispersion), d_D (dispersion slope)
    c0 = 3E8; %Speed of Light
    lambda = 1550E-9; %Pump Wavelength

    D = D*1E-12*1E9*1E-3; %Re-scale dispersion from ps/nm km to s/m m
    d2_neff = -c0*D/lambda; %Effective index of D2

    beta2 = ((lambda^3)*d2_neff/(2*pi*c0^2))*(1E12)^2; %beta2 in ps^2/m

    d_D = d_D*1E-12*1E9*1E9*1E-3; %Re-scale dispersion slope

    d3_neff = -(c0/lambda)*d_D/d2_neff; %Effective index of D3
    beta3 = -((lambda^4)/(4*pi^2*c0^3)*((lambda*d3_neff) + (3*d2_neff)))*(1E12)^3;   %beta3 in ps^3/m

end



