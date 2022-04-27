function [Pyramid_PowerMap, Pyramid_PowerMap_Align, Pyramid_PeakMap,Pyramid_PowerMap_Align_Corr] = Pyramid_v2(collisionPacket, upchirp, downchirp, SF, window_len, nfft)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明

    Power_Distribution = zeros(1, window_len*2);
    Power_Distribution(1 : window_len) = 1 : window_len;
    Power_Distribution(window_len + 1 : window_len * 2) = window_len - 1 : -1 :0;
    coeff = Power_Distribution .* [upchirp upchirp];
    coeff = conj(fliplr(coeff)); 

    Pyramid_PowerMap = zeros(2^SF, length(collisionPacket) - window_len);
    Pyramid_PowerMap_Align = zeros(2^SF, length(collisionPacket) - window_len);
    Pyramid_PeakMap = zeros(1, length(collisionPacket) - window_len);
    Pyramid_PowerMap_Align_Corr = zeros(2^SF,  length(Pyramid_PowerMap_Align(1,:)) + length(coeff) - 1);
%     Pyramid_PowerMap_Align_Corr = zeros(2^SF,  length(Pyramid_PowerMap_Align(1,:)));
    Fbin_Align = zeros(2^SF,1);
    downchirps = repmat(downchirp, 1, length(Pyramid_PowerMap_Align(1,:))/length(downchirp));
    

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
                Pyramid_PowerMap_Align_Corr(ii,:) = abs(conv(Pyramid_PowerMap_Align(ii,:),coeff));
%                 Pyramid_PowerMap_Align_Corr(ii,:) = abs(ifft(Pyramid_PowerMap_Align(ii,:).*fft(upchirp)));
%                  Pyramid_PowerMap_Align_Corr(ii,:) = (Pyramid_PowerMap_Align(ii,:).*downchirps);
    end
    
end

