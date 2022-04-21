close all;clc;
%% Load Data
filename = "frameSigs_cluster1-rx0.dat";
raw_data = read_complex_binary(filename);

%% Parameters
fs = 250e3;
SF = 8;
BW = 250e3;
% xx = raw_data(1:12000);
xx = raw_data(2095330:2107617);
SNR = 40;
xx = awgn(xx, SNR); % add gaussian noise

figure;
win = 32;
overlap = floor(win*7/8);
nfft = 2^SF*2;
[S,F,T,P] = spectrogram(xx, win, overlap, nfft, fs, 'centered');
psd = 10*log10(abs(P));

norm_freq = F/(fs/BW)/1000;     % unit: kHz
win_id = 1:1:numel(T);          % window ID

surf(win_id, norm_freq, psd);
axis tight;
shading interp;
view(0, 90);
xlim([1 numel(T)]);
ylim([-1 1]*(BW/1000)/2);
set(gca, 'YTick', (-1:0.5:1)*(BW/1000)/2);

xlabel('Win ID');
ylabel('Frequency (kHz)');
title('Spectrogram of LoRa Signal');
grid off;
set(gca, 'ticklength', [0 0]);

set(gcf, 'unit', 'centimeters', 'position', [5, 20, 40, 5]);
left_margin = 0.05;
right_margin = 0.02;
bot_margin = 0.22; 
top_margin = 0.15;
set(gca, 'position', [left_margin, bot_margin, 1-left_margin-right_margin, 1-bot_margin-top_margin]);
colorbar;
hold on;

%% Generate Upchirp and Downchirp
Ts = (2^SF)/BW;
tt = 1/fs:1/fs:Ts;
k = BW/Ts;
downchirp = exp(-1j*2*pi*(k*0.5*tt-BW/2).*tt).';
upchirp = exp(1j*2*pi*(k*0.5*tt-BW/2).*tt).';
figure;hold on;
subplot(121);
pspectrum(upchirp,fs,'spectrogram','OverlapPercent',99,'Leakage',0.85,'MinThreshold',-15,'TimeResolution',0.0001);
title("Upchirp");
subplot(122);
pspectrum(downchirp,fs,'spectrogram','OverlapPercent',99,'Leakage',0.85,'MinThreshold',-15,'TimeResolution',0.0001);
title("Downchirp");

%% LoRa PHY Decoding
status = 0;
% define status: 0,Free; 1,Preambles Detected; 2,SFD Detected
ii = 1;
last_freq = 0;
window_len = Ts * fs;
cnt = 0;
aligned_num = 1;
decode_list = [];
symbol_start_num = 1;
end_num = 1;
min_peak_power = 999999;

while ii < length(xx) - window_len
    decoding_window = xx(ii : ii + window_len - 1);
    if status == 0 % To Detect Preambles
        dechirp_signal = decoding_window .* downchirp; % get the dechirp signal
        [up_peak,freq] = max(abs(fft(dechirp_signal,2^SF)));
        if freq == last_freq
            cnt = cnt + 1;
            %%%%% To record the min power of received chirp on frequency domain
            if(cnt > 1 && abs(up_peak) < min_peak_power)
                min_peak_power = abs(up_peak);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if cnt == 4 % number of upchirps in preambles
                status = 1;
                ii = ii - freq + 1; % make decoding_window aligned with chirp signal
                aligned_num = ii; % record the aligned number
            end
        end
        last_freq = freq;

    elseif status == 1 % Preambles Detected
        dechirp_signal = decoding_window .* upchirp;
        [down_peak,~] = max(abs(fft(dechirp_signal,2^SF)));
        dechirp_signal = decoding_window .* downchirp;
        [up_peak,~] = max(abs(fft(dechirp_signal,2^SF)));
        if(abs(down_peak) > abs(up_peak)) % Downchirp Detected
            status = 2;
            previous_decoding_window = xx(ii - window_len : ii - 1);
            previous_dechirp_signal = previous_decoding_window .* upchirp;
            [previous_down_peak,~] = max(abs(fft(previous_dechirp_signal,2^SF)));
            previous_dechirp_signal = previous_decoding_window .* downchirp;
            [previous_up_peak,~] = max(abs(fft(previous_dechirp_signal,2^SF)));
            if(abs(previous_down_peak) > abs(previous_up_peak)) % Previous One is a Downchirp
                ii = ii + 1.25 * window_len;
            else
                ii = ii + 2.25 * window_len; % Start from the first symbol
            end
            symbol_start_num = ii;
            ii = ii - window_len;
        end

    elseif status == 2
        dechirp_signal = decoding_window .* downchirp;
        [peak,freq] = max(abs(fft(dechirp_signal,2^SF))); % get the dechirp signal
        if(abs(peak) > 0.9 * min_peak_power) % use the min power in preambles to judge
            decode_list = [decode_list freq];
            fprintf("winID: %d, peakPower: %.2f, bin: %d.\n", length(decode_list),abs(peak),freq)
        else
            end_num = ii; % record the end number
            status = 0;
        end
    end
    
    ii = ii + window_len;
end


%% Take a Look at Decoding Graph
ii = symbol_start_num;
decoding_graph = [];
dechirp_signal = [];
while ii < end_num - window_len
    decoding_window = xx(ii : ii + window_len - 1);
    dechirp_signal = decoding_window .* downchirp;
    temp_freq = db(abs(fft(dechirp_signal,2^SF)));
    decoding_graph = [decoding_graph temp_freq];
    ii = ii + window_len;
end

figure;
imagesc(decoding_graph);
xlabel('Symbol num');
ylabel('Decoding result');
title("Decoding Graph");
set(gca, 'YTick', (0:0.25:1)*2^SF);




