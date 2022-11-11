function [p_hat, theta] = estimate_position(x, dn_max, N, arrays)
M = size(arrays, 1);
%memory allocation
theta = zeros(M,1);
for i=1:M
    % signal of the ith array
    Y = x(:,2*i-1:2*i);
    % calculate angle of arrival
    [theta(i), ~, ~] = AOA_GCC(Y, dn_max, N);
end
% transform theta to global coordinates
theta = pi/2-(theta - [0; pi/2; pi; 3*pi/2]);
% calculate least-squares estimate of the source position
A = [tan(theta)./sqrt(tan(theta).^2+1), -1./sqrt(tan(theta).^2+1)];
b = arrays(:,2) - tan(theta).*arrays(:,1);
b = b./sqrt(tan(theta).^2+1);
p_hat = -(A'*A)\A'*b;