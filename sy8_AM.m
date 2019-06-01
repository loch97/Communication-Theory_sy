close all;clc;
% 幅度调制及解调
fm = 2000;   % 基带信号频率
fn  = 20000;  % 载波频率
N = 40001;   % 采样点数目
fs = 1e6;   % 采样率
t = 0:1/fs:(N-1)/fs;    % 时域
f = (0:(N-1))*fs/N-fs/2;    %频域
fx = cos(fm*2*pi*t);    % 基带信号
ma = 1;  % 调幅系数

% AM
figure,subplot(4,1,1);
plot(1000*t,fx);
% title('基带信号');
xlabel('时间/ms');
axis([1, 4, -1, 1]);
subplot(4,1,2);
fx_fre = (2/N)*fft(fx);
plot(f/1000,fftshift(abs(fx_fre)));  % fftshift将FFT的直流分量移到频谱中心
% title('基带频谱');
xlabel('频率/kHz');
axis([-3, 3, 0, 1]);
fs = cos(fn*2*pi*t);    %载波
subplot(4,1,3);
plot(1000*t,fs);
axis([1, 2, -1, 1]);
% title('载波信号');
xlabel('时间/ms');
subplot(4,1,4);
fs_fre = fft(fs);
plot(f/1000,(2/N)*fftshift(abs(fs_fre)));
% title('载波频谱');
xlabel('频率/kHz');
axis([-30, 30, 0, 1]);

f_am = (1+ma*fx).*fs;
figure,subplot(2,1,1);
plot(1000*t,f_am);
xlabel('调幅系数为100%的AM调制');
axis([1, 4, -2, 2]);
subplot(2,1,2);
f_am_fre = fft(f_am);
plot(f/1000,(2/N)*fftshift(abs(f_am_fre)));
xlabel('AM调制频谱');
axis([-30. 30, 0, 1])

f_am_1 = awgn(f_am, 5);
figure,subplot(4,2,1);
plot(1000*t,f_am_1);
ylabel('加噪声');title('调幅系数100%的相干解调');
axis([1, 2.5, -4, 4]);
subplot(4,2,2);
f_am_1_fre = fft(f_am_1);
plot(f/1000,(2/N)*fftshift(abs(f_am_1_fre)));
ylabel('加噪声');
axis([-40, 40, 0, 1]);
subplot(4,2,3);
fsamp = 1e6;
fcuts = [16000 17500 22500 24000];
mags = [0 1 0];
devs = [0.05 0.01 0.05];
[n, Wn, beta, ftype] = kaiserord(fcuts, mags, devs, fsamp);
hh = fir1(n, Wn, ftype, kaiser(n+1, beta), 'noscale');
f_am_bpf = fftfilt(hh, f_am_1);
plot(1000*t,f_am_bpf);
axis([1, 2.5, -2, 2]);
ylabel('经BPF');
subplot(4,2,4);
f_am_bpf_fre = fft(f_am_bpf);
plot(f/1000,(2/N)*fftshift(abs(f_am_bpf_fre)));
ylabel('经BPF');
axis([-40, 40, 0, 1]);
subplot(4,2,5);
f_am_s = 2*f_am_bpf.*fs;
plot(1000*t,f_am_s);
axis([1, 2.5, -1, 4]);
ylabel('乘载波');
subplot(4,2,6);
f_am_s_fre = fft(f_am_s);
plot(f/1000,(2/N)*fftshift(abs(f_am_s_fre)));
axis([-45, 45, 0, 2]);
ylabel('乘载波');
subplot(4,2,7);
fsamp = 1e6;
fcuts = [3000 20000];
mags = [1 0];
devs = [0.01 0.05];
[n, Wn, beta, ftype] = kaiserord(fcuts, mags, devs, fsamp);
hh = fir1(n, Wn, ftype, kaiser(n+1, beta), 'noscale');
f_am_lpf = fftfilt(hh, f_am_s);
plot(1000*t,f_am_lpf);
axis([1, 2.5, 0, 2]);
xlabel('时间/ms');
ylabel('经LPF');
subplot(4,2,8);
f_am_lpf_fre = fft(f_am_lpf);
plot(f/1000, (2/N)*fftshift(abs(f_am_lpf_fre)));
axis([-40, 40, 0, 2]);
ylabel('经LPF');
xlabel('频率/kHz');