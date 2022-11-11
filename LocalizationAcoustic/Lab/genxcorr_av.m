function [rxy, I] = genxcorr_av(x, y, N, I_old, alpha)

eps = 1e-7;
% calculate FFTs
X = fft(x(1:N));
Y = fft(y(1:N));
% calculate periodogram
I = X.*conj(Y);
% normalize periodogram
I = I./(abs(I)+eps);
% smooth
I = alpha*I_old+(1-alpha)*I;
% calculate IFFT and cut off imaginary part
rxy = real(ifft(I));
