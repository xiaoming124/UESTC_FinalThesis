close all;clc;
fs = 250e3;
SF = 8;
BW = 250e3;
%% Generate Symbol and Downchirp
Ts = (2^SF)/BW;
tt = 1/fs:1/fs:Ts;
k = BW/Ts;
windowLen = 2^SF;
nfft = 2^SF;
upchirp = exp(1j*2*pi*(k*0.5*tt-BW/2).*tt).';
downchirp = exp(-1j*2*pi*(k*0.5*tt-BW/2).*tt).';
symbol = exp(1j*2*pi*(k*0.5*tt-BW/2).*tt).';
symbol = awgn(symbol,-15); % add gaussian noise

dechirp = symbol .* downchirp;

xsignal = [zeros(windowLen, 1); symbol; zeros(windowLen, 1)];

doubleWinLen = windowLen*2;
doubleDcp = [downchirp; downchirp];
peakMap = zeros(2^SF*2,windowLen);
dechirp_fft_align = zeros(windowLen*2,1);

for ii = 1:windowLen
    dechirp_fft = fft(xsignal(ii : ii + doubleWinLen - 1) .* doubleDcp, doubleWinLen);
    idx = mod(ii,windowLen);
    for bin = 1:nfft*2
        idx = mod(nfft*2 + bin +2 - ii*2, nfft*2);
        dechirp_fft_align(1 + idx) = dechirp_fft(bin);
        peakMap(:,ii) = dechirp_fft_align';
    end
end

doubleDechirpPeakList = zeros(windowLen, 1);
for ii = 1:windowLen
    doubleDechirpPeakList(ii) = max(abs(fft(downchirp.*peakMap(ii*2,:)')));
end

figure;hold on;
subplot(211);
pspectrum(symbol,fs,'spectrogram','OverlapPercent',99,'Leakage',0.85,'MinThreshold',-15,'TimeResolution',0.0001);
title("Symbol Chirp");
subplot(212);
plot(abs(fft(dechirp)));
title("Dechirp FFT");

figure;hold on;
subplot(211);
surf(abs(peakMap));
axis tight;
shading interp;
% view(0, 0);
title('Peak Map');
grid off;
colorbar;
hold on;
subplot(212);
plot(doubleDechirpPeakList);
title('Peak Map Dechirp FFT');

