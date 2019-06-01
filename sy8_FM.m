close all;clc;
% FM调制
fm = 2000;   % 基带信号频率
fc  = 20000;  % 载波频率
N = 40001;   % 采样点数目
fs = 1e6;   % 采样率
t = 0:1/fs:(N-1)/fs;    % 时域
f = (0:(N-1))*fs/N-fs/2;    %频域
mt = cos(fm*2*pi*t);    % 基带信号/信息信号

w1 = 0; w2 = 0; %　基带信号积分
for m = 1:length(t)
    w1 = mt(m) + w2;
    w2 = mt(m) + w1;
    fi(m) = w1/(2*fs);
end
fi = fi*2*pi/max(abs(fi));
kfm = 1;A = 1;
HI = cos(kfm*fi);
HQ = sin(kfm*fi);

yo = A*cos(2*pi*fc*t).*HI - A*sin(2*pi*fc*t).*HQ;   %FM信号
figure,subplot(2,1,1);
plot(1000*t,yo);title('FM调制信号');
axis([0 2 -1 1]);
xlabel('时间/ms');ylabel('时域波形');
subplot(2,1,2);
yo_fre = 10*fft(yo);
plot(f/1000,(2/N)*fftshift(abs(yo_fre)));
xlabel('频率/kHz');ylabel('频域波形');
axis([-40 40 0 4]);

% 加噪声
yoo = awgn(yo,30);
figure,subplot(2,1,1);
plot(1000*t,yoo);title('加噪声的FM调制信号');
axis([0 2 -1 1]);
xlabel('时间/ms');ylabel('时域波形');
subplot(2,1,2);
yoo_fre = 10*fft(yoo);
plot(f/1000,(2/N)*fftshift(abs(yoo_fre)));
xlabel('频率/kHz');ylabel('频域波形');
axis([-40 40 0 4]);

%　带通滤波
fsamp = 1e6;
fcuts = [0.01 1000 40000 41000];
mags = [0 1 0];
devs = [0.05 0.01 0.05];
[n, Wn, beta, ftype] = kaiserord(fcuts, mags, devs, fsamp);
hh = fir1(n, Wn, ftype, kaiser(n+1, beta), 'noscale');
yoo_bpf = fftfilt(hh, yoo);
figure,subplot(2,1,1);
plot(1000*t,yoo_bpf);title('加噪声信号通过带通滤波器');
axis([2 6 -1 1]);
xlabel('时间/ms');ylabel('时域波形');
subplot(2,1,2);
yoo_bpf_fre = 10*fft(yoo_bpf);
plot(f/1000,(2/N)*fftshift(abs(yoo_bpf_fre)));
xlabel('频率/kHz');ylabel('频域波形');
axis([-40 40 0 4]);

% 接收信号微分处理
for i=1:length(t)-1
    diff_st_pb(i)=(yoo_bpf(i+1) - yoo_bpf(i))/(1/fs); % o_bpf o_bpf
end

sd_noise = abs(hilbert(diff_st_pb)); % 包络检波
figure,subplot(3,1,1);
plot(1000*t,[sd_noise 0]);
title('包络得到的信号');
xlabel('时间/ms');ylabel('时域波形');
axis([2, 4, 0, 0.25e6]);
fsamp = 1e6;
fcuts = [0.01 3000];
mags = [0 1];
devs = [0.01 0.05];
[n, Wn, beta, ftype] = kaiserord(fcuts, mags, devs, fsamp);
hh = fir1(n, Wn, ftype, kaiser(n+1, beta), 'noscale');
sd_out = fftfilt(hh,sd_noise/1e5);
subplot(3,1,2);
plot(1000*t, [sd_out 0]);
title('解调信号（去除直流分量）');
xlabel('时间/ms');ylabel('时域波形');
axis([5, 10, -1, 1]);
subplot(3,1,3);
sd_out_fre = fft(sd_out);
plot(f/1000,[(2/N)*fftshift(abs(sd_out_fre)) 0]);
axis([-3 3 0 1]);
xlabel('频率/kHz');ylabel('频域波形');