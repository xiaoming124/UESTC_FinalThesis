function [Pyramid_PowerMap, Pyramid_PowerMap_Align, Pyramid_PeakMap] = Pyramid(collisionPacket, downchirp, SF, window_len, nfft)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
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
end

