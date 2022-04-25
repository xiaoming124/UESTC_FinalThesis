function [Pyramid_PowerMap, Pyramid_PowerMap_Align, Pyramid_PeakMap,Pyramid_PowerMap_Align_Corr] = Pyramid(collisionPacket, upchirp, downchirp, SF, window_len, nfft)
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
    Pyramid_PowerMap = zeros(2^SF, length(collisionPacket) - window_len);
    Pyramid_PowerMap_Align = zeros(2^SF, length(collisionPacket) - window_len);
    Pyramid_PeakMap = zeros(1, length(collisionPacket) - window_len);
    Pyramid_PowerMap_Align_Corr = zeros(2^SF*2,  length(Pyramid_PowerMap_Align(1,:)) + length(upchirp) - 1);
    Fbin_Align = zeros(2^SF,1);

    for ii = 1:length(collisionPacket) - window_len
        decoding_window = collisionPacket(ii : ii + window_len - 1);
        dechirp_signal = decoding_window .* downchirp;
        Fbin = ((fft(dechirp_signal,2^SF)));
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
    
    for ii = 1 : 2^SF
                Pyramid_PowerMap_Align_Corr(ii,:) = abs(conv(Pyramid_PowerMap_Align(ii,:),fft(upchirp)));
    end
    
end

