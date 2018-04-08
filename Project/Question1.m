%% Question 1
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
xlabel('Sample'),ylabel('Signal');
legend('x','y','signal');
figure()
plot(w,abs(H));
title('Frequency Response of Adaptive Filter');
xlabel('\omega (rad/s)'),ylabel('Frequency Response');
figure()
plot(a);
title('Convergence of a');
xlabel('Sample'),ylabel('a');
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
xlabel('Sample'),ylabel('Signal');
legend('x','y','signal');
figure()
plot(w,abs(H));
title('Frequency Response of Adaptive Filter');
xlabel('\omega (rad/s)'),ylabel('Frequency Response');
figure()
plot(a);
title('Convergence of a');
xlabel('Sample'),ylabel('a');
figure()
spectrogram(y);
title('Spectra');

%%
% <latex>
% The adaptive filter is able to quickly attenuate the noise, and to track
% change in frequency.  We see this in the linear pattern in the plot of a.
% Furthermore, we can see that the moving noise band in the spectrogram plot
% is considerably attenuated, so the adaptive filter is able to find the new
% frequency and change with it.
% </latex>

%% Question 1c

w = pi/2; 
wnoise1 = 5*pi/9; 
wnoise2 = 7*pi/8;
A = 1;
Anoise = 4;
n = 0:1000;
signal = A*sin(w.*n);
noise1 = Anoise*sin(wnoise1.*n);
noise2 = Anoise*sin(wnoise2.*n);

x = signal+noise1+noise2;

r1 = 0.9;
r2 = 0.85;
mu1 = 0.0005;
mu2 = 0.001;

% First filter
e1 = zeros(1,length(x));
y1 = zeros(1,length(x));
a1 = zeros(1,length(x));

a1(1) = 0;
a1(2) = 0;
e1(1) = x(1);
y1(1) = e1(1);
e1(2) = x(2) + a1(1)*x(1);
y1(2) = e1(2) - r1*a1(1)*y1(1);

a1(3) = a1(2) - mu1*y1(2)*x(1);

for i = 3:length(n)
    e1(i) = x(i) + a1(i)*x(i-1) + x(i-2);
    y1(i) = e1(i) - r1*a1(i)*y1(i-1)-r1^2*y1(i-2);

    a1(i+1) = a1(i) - mu1*y1(i)*x(i-1);
    if a1(i+1) < -2
        a1(i+1) = 0;
    end
    if a1(i+1) > 2
        a1(i+1) = 0;
    end
end

% figure()
% stem(x,'g')
% hold on
% stem(y1,'r')
% stem(signal,'b')
% legend('x','y','signal')
% title('Notch filter adjusts to signal')
% xlabel('Sample')
% ylabel('Signal')
% hold off

% Cascaded second filter
x2 = y1;

e2 = zeros(1,length(x));
y2 = zeros(1,length(x));
a2 = zeros(1,length(x));

a2(1) = 0;
a2(2) = 0;
e2(1) = x2(1);
y2(1) = e2(1);
e2(2) = x2(2) + a2(1)*x2(1);
y2(2) = e2(2) - r2*a2(1)*y2(1);

a2(3) = a2(2) - mu2*y2(2)*x2(1);

for i = 3:length(n)
    e2(i) = x2(i) + a2(i)*x2(i-1) + x2(i-2);
    y2(i) = e2(i) - r2*a2(i)*y2(i-1)-r2^2*y2(i-2);

    a2(i+1) = a2(i) - mu2*y2(i)*x2(i-1);
    if a2(i+1) < -2
        a2(i+1) = 0;
    end
    if a2(i+1) > 2
        a2(i+1) = 0;
    end
end

figure()
stem(x,'g')
hold on
stem(y2,'r')
stem(signal,'b')
legend('x','y','signal')
title('Notch filter adjusts to signal')
xlabel('Sample')
ylabel('Signal')
hold off

w1 = linspace(-pi,pi,length(x));
X = dtft(n,w1,x);
Y = dtft(n,w1,y1);

H = Y./X;
figure()
plot(w1,H,[-7*pi/8 -7*pi/8],[min(H) max(H)],'r', ...
    [7*pi/8 7*pi/8],[min(H) max(H)], 'r')
title('Frequency response of adaptive notch filter')
xlabel('\omega (rad/s)')
ylabel('H')

figure()
spectrogram(y1)
title('High frequencies are cutoff after adaptation to noise frequency')

figure()
plot(1:length(a1),a1)
title('Convergence of coefficient a to constant value')
xlabel('Sample')
ylabel('a')

