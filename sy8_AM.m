close all;clc;
% ���ȵ��Ƽ����
fm = 2000;   % �����ź�Ƶ��
fn  = 20000;  % �ز�Ƶ��
N = 40001;   % ��������Ŀ
fs = 1e6;   % ������
t = 0:1/fs:(N-1)/fs;    % ʱ��
f = (0:(N-1))*fs/N-fs/2;    %Ƶ��
fx = cos(fm*2*pi*t);    % �����ź�
ma = 1;  % ����ϵ��

% AM
figure,subplot(4,1,1);
plot(1000*t,fx);
% title('�����ź�');
xlabel('ʱ��/ms');
axis([1, 4, -1, 1]);
subplot(4,1,2);
fx_fre = (2/N)*fft(fx);
plot(f/1000,fftshift(abs(fx_fre)));  % fftshift��FFT��ֱ�������Ƶ�Ƶ������
% title('����Ƶ��');
xlabel('Ƶ��/kHz');
axis([-3, 3, 0, 1]);
fs = cos(fn*2*pi*t);    %�ز�
subplot(4,1,3);
plot(1000*t,fs);
axis([1, 2, -1, 1]);
% title('�ز��ź�');
xlabel('ʱ��/ms');
subplot(4,1,4);
fs_fre = fft(fs);
plot(f/1000,(2/N)*fftshift(abs(fs_fre)));
% title('�ز�Ƶ��');
xlabel('Ƶ��/kHz');
axis([-30, 30, 0, 1]);

f_am = (1+ma*fx).*fs;
figure,subplot(2,1,1);
plot(1000*t,f_am);
xlabel('����ϵ��Ϊ100%��AM����');
axis([1, 4, -2, 2]);
subplot(2,1,2);
f_am_fre = fft(f_am);
plot(f/1000,(2/N)*fftshift(abs(f_am_fre)));
xlabel('AM����Ƶ��');
axis([-30. 30, 0, 1])

f_am_1 = awgn(f_am, 5);
figure,subplot(4,2,1);
plot(1000*t,f_am_1);
ylabel('������');title('����ϵ��100%����ɽ��');
axis([1, 2.5, -4, 4]);
subplot(4,2,2);
f_am_1_fre = fft(f_am_1);
plot(f/1000,(2/N)*fftshift(abs(f_am_1_fre)));
ylabel('������');
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
ylabel('��BPF');
subplot(4,2,4);
f_am_bpf_fre = fft(f_am_bpf);
plot(f/1000,(2/N)*fftshift(abs(f_am_bpf_fre)));
ylabel('��BPF');
axis([-40, 40, 0, 1]);
subplot(4,2,5);
f_am_s = 2*f_am_bpf.*fs;
plot(1000*t,f_am_s);
axis([1, 2.5, -1, 4]);
ylabel('���ز�');
subplot(4,2,6);
f_am_s_fre = fft(f_am_s);
plot(f/1000,(2/N)*fftshift(abs(f_am_s_fre)));
axis([-45, 45, 0, 2]);
ylabel('���ز�');
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
xlabel('ʱ��/ms');
ylabel('��LPF');
subplot(4,2,8);
f_am_lpf_fre = fft(f_am_lpf);
plot(f/1000, (2/N)*fftshift(abs(f_am_lpf_fre)));
axis([-40, 40, 0, 2]);
ylabel('��LPF');
xlabel('Ƶ��/kHz');