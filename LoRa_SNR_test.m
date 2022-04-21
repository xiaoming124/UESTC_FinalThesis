close all;clc;
fs = 250e3;
SF = 8;
BW = 250e3;
%% Generate Symbol and Downchirp
Ts = (2^SF)/BW;
tt = 1/fs:1/fs:Ts;
k = BW/Ts;
downchirp = exp(-1j*2*pi*(k*0.5*tt-BW/2).*tt).';
symbol = exp(1j*2*pi*(k*0.5*tt-BW/2).*tt).';
% symbol = awgn(symbol,-15); % add gaussian noise

dechirp = symbol .* downchirp;

figure;hold on;
subplot(311);
pspectrum(symbol,fs,'spectrogram','OverlapPercent',99,'Leakage',0.85,'MinThreshold',-15,'TimeResolution',0.0001);
title("Symbol Chirp under Noise");
subplot(312);
plot((abs(fft(dechirp, 2^SF))));
subplot(313);
plot(db(abs(fft(dechirp, 2^SF))));
