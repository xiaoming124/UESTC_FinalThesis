close all;clc;
fs = 250e3;
SF = 8;
BW = 250e3;
% SNR = 40;

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

%% Pyramid

 SNR = 20;
collisionPacket = [zeros(window_len,1);symbol;zeros(window_len,1)];
collisionPacket = awgn(collisionPacket, SNR);
[Pyramid_PowerMap,Pyramid_PowerMap_Align,Pyramid_PeakMap] = Pyramid(collisionPacket, downchirp, SF, window_len, nfft);

figure('Name','Pyramid_PowerMap','position',[0,500,1500,300]);

subplot(141);
surf(Pyramid_PowerMap_Align);
axis tight;
shading interp;
% view(0, 0);
title('SNR = 20dB');
grid off;
% colorbar;
hold on;

 SNR = 5;

collisionPacket = [zeros(window_len,1);symbol;zeros(window_len,1)];
collisionPacket = awgn(collisionPacket, SNR);
[Pyramid_PowerMap,Pyramid_PowerMap_Align,Pyramid_PeakMap] = Pyramid(collisionPacket, downchirp, SF, window_len, nfft);

subplot(142);
surf(Pyramid_PowerMap_Align);
axis tight;
shading interp;
% view(0, 0);
title('SNR = 5dB');
grid off;
% colorbar;
hold on;

SNR = -5;

collisionPacket = [zeros(window_len,1);symbol;zeros(window_len,1)];
collisionPacket = awgn(collisionPacket, SNR);
[Pyramid_PowerMap,Pyramid_PowerMap_Align,Pyramid_PeakMap] = Pyramid(collisionPacket, downchirp, SF, window_len, nfft);

subplot(143);
surf(Pyramid_PowerMap_Align);
axis tight;
shading interp;
% view(0, 0);
title('SNR = -5dB');
grid off;
% colorbar;
hold on;

SNR = -10;

collisionPacket = [zeros(window_len,1);symbol;zeros(window_len,1)];
collisionPacket = awgn(collisionPacket, SNR);
[Pyramid_PowerMap,Pyramid_PowerMap_Align,Pyramid_PeakMap] = Pyramid(collisionPacket, downchirp, SF, window_len, nfft);

subplot(144);
surf(Pyramid_PowerMap_Align);
axis tight;
shading interp;
% view(0, 0);
title('SNR = -10dB');
grid off;
% colorbar;
hold on;

