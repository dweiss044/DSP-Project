%% PART A
clear all;
clc;
w0=pi/2;
n=0:1:1000;
xdesired=1.*sin(w0.*n);
%stem(n,xdesired);
noise=5.*sin((7*pi/8).*n);
x=noise+xdesired;
%% part (a)
mu=0.0005;
r=0.9;
e=1:length(n);
y=1:length(n);
a=1:length(n);
a(1)=0;%-2.*cos(w0);
a(2)=a(1);
e(1)=x(1);
e(2)=x(2)+a(2).*x(1);
y(1)=e(1);
y(2)=e(2)-r.*a(2).*y(1);
a(3)=a(2)-mu.*y(2).*x(1);
for i=3:length(n)
    e(i)=x(i)+a(i).*x(i-1)+x(i-2);
    y(i)=e(i)-r.*a(i).*y(i-1)-r^2.*y(i-2);
    a(i+1)=a(i)-mu.*y(i).*x(i-1);
    if a(i+1)>2
        a(i+1)=0;
    end
    if a(i+1)<-2
         a(i+1)=0;
     end
end
w=linspace(-pi,pi,256);
X=dtft(x,w,n);
Y=dtft(y,w,n);
H=Y./X;
figure()
stem(n,x);hold on
stem(n,y);hold on
stem(n,xdesired);
title('Signal filtered and desired');
xlabel('n'),ylabel('Signal');
legend('x','y','x desired');
figure()
plot(w,abs(H));
title('Frequency Response of Adaptive Filter');
xlabel('Frequency'),ylabel('Frequency Response');
figure()
plot(a);
title('Convergence of a');
xlabel('time n'),ylabel('a');
figure()
spectrogram(y);
title('Spectra');
%% part (2)
phi=7*pi/8 + pi/10^5.*n;
xdesired=1.*sin(w0.*n);
%stem(n,xdesired);
noise=5.*sin((phi.*n));
x=noise+xdesired;
mu=0.0005;
r=0.9;
e=1:length(n);
y=1:length(n);
a=1:length(n);
a(1)=0;%-2.*cos(w0);
a(2)=a(1);
e(1)=x(1);
e(2)=x(2)+a(2).*x(1);
y(1)=e(1);
y(2)=e(2)-r.*a(2).*y(1);
a(3)=a(2)-mu.*y(2).*x(1);
for i=3:length(n)
    e(i)=x(i)+a(i).*x(i-1)+x(i-2);
    y(i)=e(i)-r.*a(i).*y(i-1)-r^2.*y(i-2);
    a(i+1)=a(i)-mu.*y(i).*x(i-1);
    if a(i+1)>2
        a(i+1)=0;
    end
    if a(i+1)<-2
         a(i+1)=0;
     end
end
w=linspace(-pi,pi,256);
X=dtft(x,w,n);
Y=dtft(y,w,n);
H=Y./X;
figure()
stem(n,x);hold on
stem(n,y);hold on
stem(n,xdesired);
title('Signal filtered and desired');
xlabel('n'),ylabel('Signal');
legend('x','y','x desired');
figure()
plot(w,abs(H));
title('Frequency Response of Adaptive Filter');
xlabel('Frequency'),ylabel('Frequency Response');
figure()
plot(a);
title('Convergence of a');
xlabel('time n'),ylabel('a');
figure()
spectrogram(y);
title('Spectra');
