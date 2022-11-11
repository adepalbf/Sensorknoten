function rxy = xycorr(x, y, N)
% calculate FFTs
X = fft(x(1:N));
Y = fft(y(1:N));
% calculate periodogram
I = X.*conj(Y)/N;
% calculate IFFT and cut off imaginary part
rxy = real(ifft(I));