function y = AddGausNoise(x, sigma)
% adds Gaussian noise with standard deviation sigma

	y = x + sigma*randn(size(x));
end