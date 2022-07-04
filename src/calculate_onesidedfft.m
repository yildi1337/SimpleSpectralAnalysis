function [f, X_onesided] = calculate_onesidedfft(x, window, NFFT, fs)
% Calculates the one-sided and complex-valued FFT of the input signal x. In 
% addition, the according frequency vector f in Hz is calculated.
%
% Note that the length of the NFFT must be an integer divisor of the length
% of the input signal x. The window needs to have the length NFFT and fs is
% the sample rate in Hz.
%
% pd (2019)

% frequency vector
f = fs .* (0 : (NFFT/2)) ./ NFFT;

% apply window
x = x .* window;

% determine two- and one-sided (complex-valued) spectra
X_twosided = fft(x, NFFT);    
X_onesided = X_twosided(1:NFFT/2+1);
