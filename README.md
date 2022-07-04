# SimpleSpectralAnalysis
A simple Matlab function for calculating various spectra/spectral densities

# Details
The function ''calculate_spectra'' (together with the subfunction ''calculate_onesidedfft'') allows to calculate the

* Power Spectral Density (PSD)
* Power Spectrum (PS)
* Linear/Amplitude Spectral Density (LSD)
* Linear/Amplitude Spectrum (LS)

for a given input (time-domain) signal. In the test script, the time-domain signal is given by a sinusoidal voltage signal with an amplitude of $1~\mathrm{V}$ and a frequency of $1~\mathrm{kHz}$ and superimposed white noise with a voltage noise density of $1~\mu\mathrm{V}/\sqrt{\mathrm{Hz}}$. **This function is working without the Signal Processing Toolbox.**

# Test Script
```matlab
% Script for testing the function calculate_spectra which calculates
% various spectra/spectral densities for a given input signal.
% In this test, the test signal x is given by a sinusoidal voltage signal
% with the frequency fx and superimposed white noise.
%
% pd (2019)
clc;
clear all;
close all;

% parameters
fs = 250e3;               % sample rate in Hz
Ts = 1/fs;                % sample period in s
Tsim = 10;                % signal length in s
fx = 1e3;                 % test signal frequency in Hz
Ax = 1;                   % amplitude of test signal in V
whitenoisex = 1e-6;       % white noise density of test signal in V/sqrt(Hz)
NFFT = fs/10;             % length of single FFT
mywindow = hann(NFFT);    % analysis window to use

% time vector
t = 0 : Ts : Tsim-Ts;

% sinusoidal test signal with noise
x = Ax .* cos(2*pi*fx*t) + sqrt(whitenoisex.^2 .* fs ./ 2) .* randn(1,length(t));

% calculate all spectra
[f, PSD, PS, LSD, LS, ENBW, Naverages] = calculate_spectra(x, mywindow, NFFT, fs);

% plot the results
figure(1);

subplot(2,2,1);
loglog(f,PSD);
xlim([10 100e3]);
ylim([1e-13 1]);
yticks(logspace(-13,0,14));
grid on;
xlabel('Frequency [Hz]');
ylabel('Power spectral density [W/Hz]');
title([ 'PSD: ENBW = ' num2str(ENBW) ' Hz, Naverages = ' num2str(Naverages) ]);

subplot(2,2,2);
loglog(f,PS);
xlim([10 100e3]);
ylim([1e-13 1]);
yticks(logspace(-13,0,14));
grid on;
xlabel('Frequency [Hz]');
ylabel('Power spectrum [W]');
title([ 'PS: ENBW = ' num2str(ENBW) ' Hz, Naverages = ' num2str(Naverages) ]);

subplot(2,2,3);
loglog(f,LSD);
xlim([10 100e3]);
ylim([1e-7 1]);
yticks(logspace(-7,0,8));
grid on;
xlabel('Frequency [Hz]');
ylabel('Amplitude spectral density [V/Hz^{0.5}]');
title([ 'LSD: ENBW = ' num2str(ENBW) ' Hz, Naverages = ' num2str(Naverages) ]);

subplot(2,2,4);
loglog(f,LS);
xlim([10 100e3]);
ylim([1e-7 1]);
yticks(logspace(-7,0,8));
grid on;
xlabel('Frequency [Hz]');
ylabel('Amplitude spectrum [V]');
title([ 'LS: ENBW = ' num2str(ENBW) ' Hz, Naverages = ' num2str(Naverages) ]);
```

# Functions
```matlab
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
```

```matlab
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
```

# Result

<p align="center">
  <img src="https://github.com/yildi1337/SimpleSpectralAnalysis/blob/main/img/results.png" />
</p>

# References

I can recommend the following articles that explain more details of spectral analysis:

* G. Heinzel, A. Rüdiger, R. Schilling, Spectrum and spectral density estimation by the Discrete Fourier transform (DFT), including a comprehensive list of window functions and some new flat-top windows, 2002. https://pure.mpg.de/pubman/faces/ViewItemOverviewPage.jsp?itemId=item_152164
* P. Welch, The use of fast Fourier transform for the estimation of power spectra: A method based on time averaging over short, modified periodograms, 1967. http://doi.org/10.1109/TAU.1967.1161901
* H. Schmid, How to use the FFT and Matlab’s pwelch function for signal and noise simulations and measurements, 2012. http://www.schmid-werren.ch/hanspeter/publications/2012fftnoise.pdf
