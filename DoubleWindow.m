close all;clc;
fs = 250e3;
SF = 8;
BW = 250e3;
SNR = 10;
%% Generate Symbol and Downchirp
Ts = (2^SF)/BW;
tt = 1/fs:1/fs:Ts;
k = BW/Ts;
window_len = Ts * fs;
nfft = 2^SF;
downchirp = exp(-1j*2*pi*(k*0.5*tt-BW/2).*tt);
upchirp = exp(1j*2*pi*(k*0.5*tt-BW/2).*tt);

symbol1 = [exp(1j*2*pi*(k*0.5*tt-BW*3/8).*tt).' ; zeros(window_len,1)];
symbol2 = [zeros(window_len,1) ; exp(1j*2*pi*(k*0.5*tt+BW/4).*tt).'];
symbol3 = [zeros(window_len/2,1) ; exp(1j*2*pi*(k*0.5*tt+BW/4).*tt).' ; zeros(window_len/2,1)];
symbol = symbol1 + symbol2 + symbol3;
% symbol = symbol1;

%Peak Location: (33,257) (129,512) (193,385)

%% Pyramid
collisionPacket = [zeros(window_len,1);symbol;zeros(window_len,1)].';
collisionPacket = awgn(collisionPacket, SNR);
[Pyramid_PowerMap,Pyramid_PowerMap_Align,Pyramid_PeakMap, Pyramid_PowerMap_Align_Corr] = Pyramid_v2(collisionPacket, upchirp, downchirp, SF, window_len, nfft);

figure('position',[500,500,500,500]);
subplot(311);
surf(abs(Pyramid_PowerMap));
axis tight;
shading interp;
% view(0, 0);
title('Pyramid Power Map');
grid off;
colorbar;
hold on;

subplot(312);
surf(abs(Pyramid_PowerMap_Align));
axis tight;
shading interp;
% view(0, 0);
title('Pyramid Power Map - Aligned');
grid off;
colorbar;
hold on;

subplot(313);
surf((Pyramid_PowerMap_Align_Corr));
axis tight;
shading interp;
% view(0, 0);
title('Pyramid Power Map - Accumulated');
grid off;
colorbar;
hold on;


% disp("Pyramid Decode Result");
% disp("Ground Truth:Peak Location: (33,257) (129,512) (193,385)");
% [~, time] = max(Pyramid_PowerMap_Align(33,:));
% disp(["Aligned_Freq" 33 "Peak_Time" time]);
% [~, time] = max(Pyramid_PowerMap_Align(129,:));
% disp(["Aligned_Freq" 129 "Peak_Time" time]);
% [~, time] = max(Pyramid_PowerMap_Align(193,:));
% disp(["Aligned_Freq" 193 "Peak_Time" time]);

%% DoubleWindow
% Power_Distribution = zeros(1, window_len);
% Power_Distribution(1 : window_len) = 1 : window_len;
% Power_Distribution(window_len : window_len * 2) = window_len;
% Power_Distribution(window_len * 2 + 1 : window_len * 3) = window_len : -1 :1;

tt = 1/fs:1/fs:2*Ts;
k = BW/Ts;
Double_window_len = Ts * fs * 2;
doubleDownchirp = exp(-1j*2*pi*(k*0.5*tt-BW/2).*tt);

collisionPacket = [zeros(Double_window_len,1);symbol;zeros(Double_window_len,1)].';
collisionPacket = awgn(collisionPacket, SNR);
[DW_PowerMap, DW_PowerMap_Align, DW_PeakMap,DW_PowerMap_Align_Corr] = DoubleWin_v2(collisionPacket, upchirp, doubleDownchirp, SF, Double_window_len, nfft);

% DW_PowerMap_Align_Accumulate = zeros(2^SF*2, length(collisionPacket) - Double_window_len - 2^SF);
% DW_PowerMap_Align_Corr = zeros(2^SF*2, 2 * length(DW_PowerMap_Align(1,:)) - 1);

% for ii = 1 : 2^SF*2
%     for jj = 1 : length(DW_PowerMap_Align) - 2^SF
%         DW_PowerMap_Align_Accumulate(ii,jj) = sum(DW_PowerMap_Align(ii,jj:jj+2^SF));
%     end
% end

% for ii = 1 : 2^SF*2
%         DW_PowerMap_Align_Corr(ii,:) = xcorr(DW_PowerMap_Align(ii,:), Power_Distribution)(:,768:end);
% end

% figure('position',[1000,500,500,500]);
% subplot(311);
% surf(abs(DW_PowerMap));
% axis tight;
% shading interp;
% % view(0, 0);
% title('DoubleWindow Power Map');
% grid off;
% colorbar;
% hold on;
% 
% subplot(312);
% surf(abs(DW_PowerMap_Align));
% axis tight;
% shading interp;
% % view(0, 0);
% title('DoubleWindow Power Map - Aligned');
% grid off;
% colorbar;
% hold on;
% 
% 
% subplot(313);
% % surf(DW_PowerMap_Align_Accumulate);
% surf(abs(DW_PowerMap_Align_Corr));
% axis tight;
% shading interp;
% % view(0, 0);
% title('DoubleWindow Power Map - Accumulated');
% grid off;
% colorbar;
% hold on;

% disp("DoubleWindow Decode Result");
% disp("Ground Truth:Peak Location: (64,257) (256,512) (383,385)");
% [~, time] = max(DW_PowerMap_Align_Accumulate(64,:));
% disp(["Aligned_Freq" 64 "Peak_Time" time]);
% [~, time] = max(DW_PowerMap_Align_Accumulate(256,:));
% disp(["Aligned_Freq" 256 "Peak_Time" time]);
% [~, time] = max(DW_PowerMap_Align_Accumulate(383,:));
% disp(["Aligned_Freq" 384 "Peak_Time" time]);

% 
% disp("DoubleWindow Decode Result");
% disp("Ground Truth:Peak Location: (66,257) (258,512) (386,385)");
% [~, time] = max(DW_PowerMap_Align_Corr(66,:));
% disp(["Aligned_Freq" 66 "Peak_Time" time]);
% [~, time] = max(DW_PowerMap_Align_Corr(258,:));
% disp(["Aligned_Freq" 258 "Peak_Time" time]);
% [~, time] = max(DW_PowerMap_Align_Corr(386,:));
% disp(["Aligned_Freq" 386 "Peak_Time" time]);


% xx = DW_PowerMap_Align_Corr(130, 518:518+255);
%  xx = Pyramid_PowerMap_Align_Corr(33, :);
xx = Pyramid_PowerMap_Align(33, :);
% xx = DW_PowerMap_Align_Corr(66, :);

% FinalMap = zeros(1, length(xx));
% FinalMap2 = zeros(1, length(xx));
% for ii = 1:length(xx) - window_len
%     [peak , loc] = max(abs(fft(xx(ii:ii+window_len))));
%     FinalMap(ii) = peak;
%     FinalMap2(ii) = loc;
% end
% % xx = Pyramid_PowerMap_Align_Corr(33, 1:512);
% figure;
% plot(FinalMap);
% figure;
% plot(FinalMap2);
% [~, time] = max(FinalMap);
% disp([ time]);

% figure;
% plot(abs(xx))
tt = 1/fs:1/fs:Ts;
chirpList = zeros(1, 256);
num = 0;
symbol4 = [zeros(1,256) exp(1j*2*pi*(k*0.5*tt-BW/2+BW*num/256).*tt) zeros(1,256)];
for ii = 1:512
    yy = symbol4(ii:ii+255) .* downchirp;
    tmp = fft(yy);
    chirpList(ii) = tmp(mod(ii+num-1,256)+1);
end
figure;
subplot(211);
plot(abs(chirpList));
subplot(212);
plot(phase(chirpList));

% figure;
% plot(phase(symbol4(1:256)));

% figure;hold on;
% plot(phase(exp(1j*2*pi*(k*0.5*tt-BW/2+BW*num/256).*tt)));
% plot(phase(fft(exp(1j*2*pi*(k*0.5*tt-BW/2+BW*num/256).*tt))));
% legend('time','freq');


% zz = exp(1j*2*pi*(k*0.5*tt-BW/2+BW*num/256).*tt) .* downchirp;
% tt = 0:255;
% zz = exp(1j*pi/4*tt);
% zz2 = angle(fft(zz));
% figure;
% subplot(211);
% plot(abs(fft(zz)));
% subplot(212);
% plot(angle(fft(zz)));
% 
% disp(zz2(33));


