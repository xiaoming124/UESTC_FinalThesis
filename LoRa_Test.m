clear;clc;close all;
filename = 'rftest';%TX_sample
count=10e5;
vec=read_complex_binary(filename, count);
% plotdat(vec,0);

%take the real part of the complex sequence
vec_imag=real(vec);

%lora packets segmentation
startt=[];
endt=[];
flag=0;
upfac=1;
s=3e-3; 
e=1e-3; 
for i=upfac*300+1:length(vec_imag)-upfac*300
    if (max(abs(vec_imag(i:i+upfac*300)))>s)&&(max(abs(vec_imag(i-upfac*300:i)))<e)&&(flag==0)
%         packet=a(i:i+6400);
%         i=i+7200;
        startt=[startt,i];
        flag=1;
    elseif ((max(abs(vec_imag(i:i+upfac*300)))<e)&&flag==1)%&&(sum(abs(a(i:i-50)))/50>5e-2)
        endt=[endt,i];
        flag=0;
    end
end


%take a look at a lora packet
fs = 500e3;
chirp_dur=480;%99
x1=vec(startt(1):endt(1));

SF = 7;
BW = 125e3;
Ts = (2^SF)/BW;
k = BW/Ts;
tt = 0:1/fs:Ts;
downchirp = exp(-1j*2*pi*(k*0.5*tt-62500).*tt).';
upchirp = exp(1j*2*pi*(k*0.5*tt-62500).*tt).';

figure;
subplot(311);

pspectrum(x1,fs,'spectrogram', ...
        'OverlapPercent',50,'Leakage',0.85,'MinThreshold',-35,'TimeResolution',.00006);

reference = [];
for ii = 1:31
    reference = [reference; upchirp];
end

x1_pad = zeros(395,1);
x1 = [x1(300:end);x1_pad];


plus_positive = x1.*reference;
subplot(312);
pspectrum(plus_positive,fs,'spectrogram', ...
        'OverlapPercent',50,'Leakage',0.85,'MinThreshold',-35,'TimeResolution',.00006);

reference = [];
for ii = 1:31
    reference = [reference; downchirp];
end


plus_negative = x1.*reference;
subplot(313);
pspectrum(plus_negative,fs,'spectrogram', ...
        'OverlapPercent',50,'Leakage',0.85,'MinThreshold',-35,'TimeResolution',.00006);

