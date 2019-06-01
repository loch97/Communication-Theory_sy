Fs=1e4;                        %����Ƶ��
len=20;                        %��Ԫ����
in=randint(1,len,4);         %������ʼ��Ԫ����
sig=[];	
out=[];
for t=1:2000              %���������ź�
    n=fix(t/100);
    if n==0
        in_a(t)=0;
    else
        in_a(t)=in(n);
    end
end
subplot(2,1,1);              %�����ź�
plot(in_a,'LineWidth',3);
title('�����ź�','FontWeight','bold','FontSize',20);
xlabel('t/s','FontSize',18);
axis([100,2100,-0.5,3.5]);
set(gca,'XTick',0:100:2000);
grid on;
cxn=xcorr(in_a,'unbiased'); %�������е�����غ���
nfft=1024;
CXk=fft(cxn,nfft);
Pxx=abs(CXk);
index=0:round(nfft/2-1);
k=index*Fs/nfft;
subplot(2,1,2);
plot_Pxx=10*log10(Pxx(index+1));
plot(k,plot_Pxx,'LineWidth',2);
title('�����źŹ�����','FontWeight','bold','FontSize',20);
axis([0,5000,-10,40]);
xlabel('Hz','FontSize',18,'FontSize',18);
for i=1:len        %���������Թ������ź�
    if  in(i)==0
        ins=[0,0];
    elseif in(i)==1
        ins=[1,0];
    elseif in(i)==2
        ins=[2,0];
    elseif in(i)==3
        ins=[3,0];
    end
    sig=[sig,ins];
end
for t=1:4000
    n=fix(t/100);
    if n==0
        s(t)=0;
    else
        s(t)=sig(n);
    end
end
figure;
subplot(2,1,1);            %�����Թ�����
plot(s,'LineWidth',3);
title('�����Թ�����','FontWeight','bold','FontSize',20);
xlabel('t/s','FontSize',18);
axis([100,4100,-0.5,3.5]);
set(gca,'XTick',0:200:4100);
grid on;
cxn=xcorr(s,'unbiased');
nfft=1024;
CXk=fft(cxn,nfft);
Pxx=abs(CXk);
index=0:round(nfft/2-1);
k=index*Fs/nfft;
subplot(2,1,2);
plot_Pxx=10*log10(Pxx(index+1));
plot(k,plot_Pxx,'LineWidth',2);
title('�����Թ����빦����','FontWeight','bold','FontSize',20);
axis([0,5000,-10,40]);
xlabel('Hz','FontSize',18);
s1=awgn(s,20);              %�������
figure;
subplot(2,1,1);
plot(s1);
title('�����������ź�','FontWeight','bold','FontSize',20);
xlabel('t/s','FontSize',18);
axis([100,4100,-0.5,3.5]);
set(gca,'XTick',0:500:4100);
cxn=xcorr(s1,'unbiased');
nfft=1024;
CXk=fft(cxn,nfft);
Pxx=abs(CXk);
index=0:round(nfft/2-1);
k=index*Fs/nfft;
subplot(2,1,2);
plot_Pxx=10*log10(Pxx(index+1));
plot(k,plot_Pxx,'LineWidth',2);
title('�����������źŹ�����','FontWeight','bold','FontSize',20);
axis([0,5000,-10,40]);
xlabel('Hz','FontSize',18);    %�˲������
fp=500;                            %ͨ����ֹ
fs= 550;                        %�����ֹ
ws=fs*2/Fs;
wp=fp*2/Fs;
[N, Wp] = ellipord(wp,ws,1,40);
[b,a]=ellip(N,1,40,Wp);
sf0=filter(b,a,s1) ;    %�˵�������������ź�
figure;
subplot(2,1,1);
plot(sf0);
axis([100,4100,-0.5,3.5]);
title('�˵�������������ź�','FontWeight','bold','FontSize',20);
xlabel('t/s','FontSize',18);
set(gca,'XTick',0:500:4100);
cxn=xcorr(sf0,'unbiased'); nfft=1024;
CXk=fft(cxn,nfft);
Pxx=abs(CXk);
index=0:round(nfft/2-1);
k=index*Fs/nfft;
subplot(2,1,2);
plot_Pxx=10*log10(Pxx(index+1));
plot(k,plot_Pxx,'LineWidth',2);
title('�˵�������������źŹ�����','FontWeight','bold','FontSize',20);
axis([0,5000,-10,40]);
xlabel('Hz','FontSize',18,'FontSize',18);
for m=1:20                  %�����о�
    p=round(sf0(200*m-20));
    out=[out,p];
end
figure;
subplot(2,1,1);        %ԭʼ��Ԫ����
stairs(in,'LineWidth',3);
title('ԭʼ��Ԫ����','FontWeight','bold','FontSize',20);
axis([1,21,-0.5,3.5]);
set(gca,'XTick',0:1:20);
grid on;
subplot(2,1,2);  %�����о���ָ����ź�����
stairs(out,'LineWidth',3);
title('�����о���ָ����ź�����','FontWeight','bold','FontSize',20);
axis([1,21,-0.5,3.5]);
set(gca,'XTick',0:1:20);
grid on;
