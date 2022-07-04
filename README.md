# SimpleSpectralAnalysis
A simple Matlab function for calculating various spectra/spectral densities

# Details
The function ''calculate_spectra'' (together with the subfunction ''calculate_onesidedfft'') allows to calculate the
* Power Spectral Density (PSD)
* Power Spectrum (PS)
* Linear/Amplitude Spectral Density (LSD)
* Linear/Amplitude Spectrum (LS)
for a given input (time-domain) signal. In the test script, the time-domain signal x is given by a sinusoidal voltage signal with a frequency of $1 \mathrm{kHz}$ and superimposed white noise with a voltage noise density of <img src="https://render.githubusercontent.com/render/math?math=1~\mathrm{V}/\sqrt{\mathrm{Hz}}">.

# Pictures

<p align="center">
  <img src="https://github.com/yildi1337/SimpleSpectralAnalysis/blob/main/img/results.png" />
</p>
