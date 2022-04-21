close all;clc;
fs = 250e3;
SF = 8;
BW = 250e3;
SNR = -5;
%% Generate Symbol and Downchirp
Ts = (2^SF)/BW;
tt = 1/fs:1/fs:Ts;
k = BW/Ts;
window_len = Ts * fs;
nfft = 2^SF;
downchirp = exp(-1j*2*pi*(k*0.5*tt-BW/2).*tt).';
% upchirp = exp(1j*2*pi*(k*0.5*tt-BW/2).*tt).';

symbol1 = [exp(1j*2*pi*(k*0.5*tt-BW*3/8).*tt).' ; zeros(window_len,1)];
symbol2 = [zeros(window_len,1) ; exp(1j*2*pi*(k*0.5*tt).*tt).'];
symbol3 = [zeros(window_len/2,1) ; exp(1j*2*pi*(k*0.5*tt-BW/4).*tt).' ; zeros(window_len/2,1)];
symbol = symbol1 + symbol2 + symbol3;
% symbol = awgn(symbol,-15); % add gaussian noise

% figure('position',[0,500,500,400]);
% subplot(211);
% pspectrum(symbol,fs,'spectrogram','OverlapPercent',99,'Leakage',0.85,'MinThreshold',-15,'TimeResolution',0.0001);
% title("Symbol Chirp under Noise");

%% Pyramid
collisionPacket = [zeros(window_len,1);symbol;zeros(window_len,1)];
collisionPacket = awgn(collisionPacket, SNR);
Pyramid_PowerMap = zeros(2^SF, length(collisionPacket) - window_len);
Pyramid_PowerMap_Align = zeros(2^SF, length(collisionPacket) - window_len);
Pyramid_PeakMap = zeros(1, length(collisionPacket) - window_len);
Fbin_Align = zeros(2^SF,1);

for ii = 1:length(collisionPacket) - window_len
    decoding_window = collisionPacket(ii : ii + window_len - 1);
    dechirp_signal = decoding_window .* downchirp;
    Fbin = (abs(fft(dechirp_signal,2^SF)));
    idx = mod(ii,window_len);
    freq_shift = mod(ii + window_len,window_len);
    for bin = 1:nfft
        idx = mod(nfft + bin - freq_shift, nfft);
        Fbin_Align(idx+1) = Fbin(bin);
    end
    Pyramid_PowerMap(:,ii) = Fbin;
    Pyramid_PowerMap_Align(:,ii) = Fbin_Align;
    Pyramid_PeakMap(:,ii) = max(Fbin);
end

% subplot(212);
% plot(Pyramid_PeakMap);
% title('Pyramid Peak Map');

figure('position',[500,500,500,600]);
subplot(211);
surf(Pyramid_PowerMap);
axis tight;
shading interp;
% view(0, 0);
title('Pyramid Power Map');
grid off;
colorbar;
hold on;

subplot(212);
surf(Pyramid_PowerMap_Align);
axis tight;
shading interp;
% view(0, 0);
title('Pyramid Power Map - Aligned');
grid off;
colorbar;
hold on;
% [Peak, Symb] = max(Pyramid_PowerMap(:,386))

%% DoubleWindow
tt = 1/fs:1/fs:2*Ts;
k = BW/Ts;
window_len = Ts * fs * 2;
doubleDownchirp = exp(-1j*2*pi*(k*0.5*tt-BW/2).*tt).';

collisionPacket = [zeros(window_len,1);symbol;zeros(window_len,1)];
collisionPacket = awgn(collisionPacket, SNR);

DW_PowerMap = zeros(2^SF*2, length(collisionPacket) - window_len);
DW_PowerMap_Align = zeros(2^SF*2, length(collisionPacket) - window_len);
DW_PeakMap = zeros(1, length(collisionPacket) - window_len);
% Fbin_Align = zeros(2^SF*2,1);

for ii = 1:length(collisionPacket) - window_len
    decoding_window = collisionPacket(ii : ii + window_len - 1);
    dechirp_signal = decoding_window .* doubleDownchirp;
    Fbin = (abs(fft(dechirp_signal,2^SF * 2)));
    idx = mod(ii,window_len);
    freq_shift = mod(ii + window_len,window_len);
    for bin = 1:nfft*2
        idx = mod(nfft*2 + bin - freq_shift*2, nfft*2);
        Fbin_Align(idx+1) = Fbin(bin);
    end
    DW_PowerMap(:,ii) = Fbin;
    DW_PowerMap_Align(:,ii) = Fbin_Align;
    DW_PeakMap(:,ii) = max(Fbin);
end

figure('position',[1000,500,500,600]);
subplot(211);
surf(DW_PowerMap);
axis tight;
shading interp;
% view(0, 0);
title('DoubleWindow Power Map');
grid off;
colorbar;
hold on;

subplot(212);
surf(DW_PowerMap_Align);
axis tight;
shading interp;
% view(0, 0);
title('DoubleWindow Power Map - Aligned');
grid off;
colorbar;
hold on;



%% Plot original signal & Peak Map
figure('position',[0,500,500,500]);
subplot(311);
pspectrum(symbol,fs,'spectrogram','OverlapPercent',99,'Leakage',0.85,'MinThreshold',-15,'TimeResolution',0.0001);
title("Symbol Chirp under Noise");

subplot(312);
plot(Pyramid_PeakMap);
title('Pyramid Peak Map');

subplot(313);
plot(DW_PeakMap);
title('DoubleWindow Peak Map');