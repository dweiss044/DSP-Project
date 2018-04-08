function [H] = dtft(n,omega,r) 
H = exp(-1i*omega'*n)*r';

end