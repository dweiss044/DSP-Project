%% Final project

clear; clc;

%% Question 1a

w = pi/2; % Hz
wnoise1 = 7*pi/8; % Hz
A = 1;
Anoise = 4;
n = [0:250];
signal = A*sin(w.*n);
noise1 = Anoise*sin(wnoise1.*n);

x = signal+noise1;

r = 0.90;
mu = 0.005;
e1 = zeros(1,length(x));
y1 = zeros(1,length(x));
a1 = zeros(1,length(x));

a1(1) = 0;
a1(2) = 0;
e1(1) = x(1);
y1(1) = e1(1);
e1(2) = x(2) + a1(1)*x(1);
y1(2) = e1(2) - r*a1(1)*y1(1);

a1(3) = a1(2) - mu*y1(2)*x(1);

for i = 3:length(n)
    e1(i) = x(i) + a1(i)*x(i-1) + x(i-2);
    y1(i) = e1(i) - r*a1(i)*y1(i-1)-r^2*y1(i-2);

    a1(i+1) = a1(i) - mu*y1(i)*x(i-1);
    if a1(i+1) < -2
        a1(i+1) = 0;
    end
    if a1(i+1) > 2
        a1(i+1) = 0;
    end
end


figure()
stem(x,'g')
hold on
stem(y1,'r')
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


%% Question 1b

w = pi/4; % Hz
wnoise1 = @(n0) 3*pi/4 + pi/1e6*n0; % Hz
A = 2;
Anoise = 5;
n = [0:2000];
signal = A*sin(w.*n);
noise1 = Anoise*sin(wnoise1(n).*n);

x = signal+noise1;

r = 0.90;
mu = 0.005;
e1 = zeros(1,length(x));
y1 = zeros(1,length(x));
a1 = zeros(1,length(x));

a1(1) = -2*cos(w);
a1(2) = a1(1);
e1(1) = x(1);
y1(1) = e1(1);
e1(2) = x(2) + a1(1)*x(1);
y1(2) = e1(2) - r*a1(1)*y1(1);

a1(3) = a1(2) - mu*y1(2)*x(1);

for i = 3:length(n)
    e1(i) = x(i) + a1(i)*x(i-1) + x(i-2);
    y1(i) = e1(i) - r*a1(i)*y1(i-1)-r^2*y1(i-2);

    a1(i+1) = a1(i) - mu*y1(i)*x(i-1);
    if a1(i+1) < -2
        a1(i+1) = 0;
    end
    if a1(i+1) > 2
        a1(i+1) = 0;
    end
end


figure()
stem(x,'g')
hold on
stem(y1,'r')
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
plot(w1,H)
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

%%
% <latex>
% The adaptive filter is able to quickly attenuate the noise, but is slow
% in tracking the change in frequency.  We see this in the oscillations of
% the convergence plot.  Because the frequency content in x is changing
% constanly, the filter cannot fully attenuate the noise, but does
% significantly reduce it.
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

