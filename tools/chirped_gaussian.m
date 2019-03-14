function u = chirped_gaussian(t, t0, FWHM, P0, w0, C)
%Chirped Gaussian Pulse
% Equation taken from: https://www.brown.edu/research/labs/mittleman/sites/brown.edu.research.labs.mittleman/files/uploads/lecture6_0.pdf
% t = time vector
% t0 = time centroid
% FWHM = intensity full width at half maximum 
% P0 = Maximum power over domain [0 1]
% w0 = carrier frequency of modulation
% C = chirp parameter

%u = P0.*exp(-1.38.*((t-t0)./FWHM).^2).*exp(1i.*(w0.*t + (C.*t.^2)));
u = P0.*exp(-1.38.*((t-t0)./FWHM).^2).*exp(1i.*(w0.*t + (C.*t.^2)));
end

