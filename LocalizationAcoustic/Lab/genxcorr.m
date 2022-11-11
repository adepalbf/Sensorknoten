function [rxy, I] = genxcorr(x, y, N)

eps = 1e-7;
% calculate FFTs
X = fft(x, N);
Y = fft(y, N);
% calculate periodogram
I = X.*conj(Y);
% normalize periodogram
I = I./(abs(I)+eps);
% calculate IFFT and cut off imaginary part
rxy = real(ifft(I));
