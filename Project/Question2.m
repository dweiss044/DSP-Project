%% Question 2

clear; clc;

%% Part 2a

for i = 1:1000
    tmp = rand;
    if tmp < 0.5
        s(i) = -1;
    else
        s(i) = 1;
    end
end

h = [0.3 1 0.7 0.3 0.2];
SNR = 25; % dB
sigma = 1/(10^(SNR/20));
w = sigma*randn(1,length(s));
x = conv(s,h,'same') + w;


