function [rxy, I] = genxcorr_av(x, N, I_old, alpha)
M = size(x,2);
eps = 1e-7;

% calculate FFTs
X = fft(x, N);
I = zeros(N,M-1);
for ii = 1:M-1
    % calculate periodogram
    I(:,ii) = conj(X(:,1)).*X(:,ii+1);
end

% normalize periodogram
I = I./(abs(I)+eps);
% smooth
I = alpha*I_old+(1-alpha)*I;
% calculate IFFT and cut off imaginary part
rxy = real(ifft(I));