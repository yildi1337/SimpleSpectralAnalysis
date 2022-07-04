function [f, PSD, PS, LSD, LS, ENBW, Naverages] = calculate_spectra(x, window, NFFT, fs)
% Calculates the
%   - Power Spectral Density (PSD)
%   - Power Spectrum (PS)
%   - Linear/Amplitude Spectral Density (LSD)
%   - Linear/Amplitude Spectrum (LS)
% for the given time domain input signal x. In addition, the according 
% frequency vector in Hz  f, the equivalent noise bandwidth (ENBW), and the
% number of averages are calculated.
%
% Note that the length of the NFFT must be an integer divisor of the length
% of the input signal x. The window needs to have the length NFFT and fs is
% the sample rate in Hz.
%
% pd (2019)

% calculate number of averages
Naverages = length(x) / NFFT;

% check if NFFT fits perfectly into the length of the input vector x
if Naverages ~= floor(Naverages)
    disp('NFFT must be an integer divisor of the length of x.');
    return
end

% split input signal into matrix (each column represents one time period)
x = reshape(x, [NFFT Naverages]);

% calculate window sums and ENBW
S1 = sum(window);
S2 = sum(window.^2);
ENBW = fs * S2/(S1^2);

% final (averaged) power spectral density
PSD = 0;

% calculate spectrum for each time period
for n = 1 : Naverages
  
    % calculate one-sided power spectrum of this time period
    [f, Xn_onesided] = calculate_onesidedfft(x(:,n), window, NFFT, fs);
    PSn = 2 * abs(Xn_onesided).^2 ./ (S1.^2);
    
    % calculate one-sided power spectral density of this time period
    PSDn = PSn ./ ENBW;
   
    % add power spectral density of this time period to the overall result
    PSD = PSD + PSDn;

end

% calculate overall (averaged) power spectral density
PSD = PSD ./ Naverages;

% calculate overall (averaged) power spectrum
PS = PSD .* ENBW;

% calculate overall (averaged) linear spectral density
LSD = sqrt(PSD);

% calculate overall (averaged) linear spectrum
LS = sqrt(PS);

