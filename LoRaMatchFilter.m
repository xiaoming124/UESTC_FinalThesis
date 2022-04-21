close all;clc;
%% Parameters
fs = 250e3;
SF = 8;
BW = 250e3;
fft_n = 2^SF;

%% Generate Upchirp and Downchirp
Ts = (2^SF)/BW;
tt = 1/fs:1/fs:Ts;
k = BW/Ts;
symbol1 = 128;
symbol2 = 250;
downchirp = exp(-1j*2*pi*(k*0.5*tt-BW/2).*tt).';
upchirp = exp(1j*2*pi*(k*0.5*tt-BW/2).*tt).';

sig1 = upchirp.*(exp(1j*2*pi*symbol1*BW/(2^SF).*tt).');
sig2 = upchirp.*(exp(1j*2*pi*symbol2*BW/(2^SF).*tt).');
sig = sig1 + sig2;
SNR = 20;
sig = awgn(sig, SNR);
%% Match Filter
coeff = conj(fliplr(upchirp)); %翻转共轭

pc_res = ifft(fft(sig,fft_n).*fft(coeff,fft_n)); % 未截取不完全滤波点
pc_res = flip(pc_res);
figure;
subplot(311);
pspectrum(sig,fs,'spectrogram','OverlapPercent',99,'Leakage',0.85,'MinThreshold',-15,'TimeResolution',0.0001);
title("Upchirp");

subplot(312);
plot(db(abs(pc_res)/max(abs(pc_res))),'r');  title('Match Filter');

subplot(313);
dcp = sig .* downchirp;
plot(db(abs(fft(dcp,fft_n))/max(abs(fft(dcp,fft_n)))));
title("LoRa Decoder");

