function [DW_PowerMap, DW_PowerMap_Align, DW_PeakMap,DW_PowerMap_Align_Corr] = DoubleWin_v2(collisionPacket, upchirp, doubleDownchirp, SF, double_window_len, nfft)
%DOUBLEWIN 此处显示有关此函数的摘要
%   此处显示详细说明
    window_len = 2^SF;
    Power_Distribution = zeros(1, window_len);
    Power_Distribution(1 : window_len) = 1 : window_len;
    Power_Distribution(window_len : window_len * 2) = window_len;
    Power_Distribution(window_len * 2 + 1 : window_len * 3) = window_len : -1 :1;
    
    DW_PowerMap = zeros(2^SF*2, length(collisionPacket) - double_window_len);
    DW_PowerMap_Align = zeros(2^SF*2, length(collisionPacket) - double_window_len);
    DW_PeakMap = zeros(1, length(collisionPacket) - double_window_len);
%     DW_PowerMap_Align_Corr = zeros(2^SF*2,  length(DW_PowerMap_Align(1,:)) + length(upchirp) - 1);
    DW_PowerMap_Align_Corr = zeros(2^SF*2,  length(DW_PowerMap_Align(1,:)));
    downchirps = repmat(doubleDownchirp, 1, length(DW_PowerMap_Align(1,:))/length(doubleDownchirp));
    % Fbin_Align = zeros(2^SF*2,1);

    for ii = 1:length(collisionPacket) - double_window_len
        decoding_window = collisionPacket(ii : ii + double_window_len - 1);
        dechirp_signal = decoding_window .* doubleDownchirp;
        Fbin = ((fft(dechirp_signal,2^SF * 2)));
        idx = mod(ii,double_window_len);
        freq_shift = mod(ii + double_window_len,double_window_len);
%         for bin = 1:nfft
%             idx = mod(nfft*2 + bin*2 - freq_shift*2, nfft*2);
%             Fbin_Align(idx+2) = Fbin(bin*2);
%         end
        
        for bin = 1:nfft*2
            idx = mod(nfft*2 + bin + 2 - freq_shift*2, nfft*2);
            Fbin_Align(idx+1) = Fbin(bin);
        end
        
        DW_PowerMap(:,ii) = (Fbin);
        DW_PowerMap_Align(:,ii) = (Fbin_Align);
        DW_PeakMap(:,ii) = max(Fbin);
    end
    
    for ii = 1 : 2^SF*2
%             tmp = xcorr(DW_PowerMap_Align(ii,:), Power_Distribution);
%             DW_PowerMap_Align_Corr(ii,:) = tmp;
%                 DW_PowerMap_Align_Corr(ii,:) = abs(conv(DW_PowerMap_Align(ii,:),fft(upchirp)));
                DW_PowerMap_Align_Corr(ii,:) = (DW_PowerMap_Align(ii,:).*downchirps);

    end

end

