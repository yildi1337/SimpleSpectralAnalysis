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

