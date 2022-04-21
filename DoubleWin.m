function [DW_PowerMap, DW_PowerMap_Align, DW_PeakMap] = DoubleWin(collisionPacket, doubleDownchirp, SF, window_len, nfft)
%DOUBLEWIN 此处显示有关此函数的摘要
%   此处显示详细说明
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
end

