function alpha = atten_to_alpha(atten)

% Convert fiber attenuation to loss coefficient alpha
%
% INPUT
%
% atten - fiber attenuation, dB/km
%
% OUTPUT
%
% alpha - loss coefficient, 1/m

alpha = atten * 1E-3 / 4.34;

end
