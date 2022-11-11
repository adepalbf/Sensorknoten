% capture time in seconds
T = 10;
% sampling frequency in hertz
f=48000;
% size of the data
N = f*T;
% record data
x=nktp_rec(N, f);
% microphone distance in meter
d = 0.175;
% speed of sound in meter/second
v = 343;
% maximum TDOA in samples
dn_max = round(f*d/v);
% array positions in meter
arrays = [2.69, 0.19];%; 5.19, 3.21; 2.70, 5.27; 0.37, 3.54];
% number of estimates to calculate, should be an integer fraction of N
M = 100;
% memory allocation
p_hat = zeros(2,M);
for t = 1:M
    x_t = x((t-1)*N/M+1:t*N/M,:);
    p_hat(:,t)  = estimate_position(x_t, dn_max, N, arrays);
end

%scatter plot to visualize movement
scatter(p_hat(1,:), p_hat(2,:))
hold on
scatter(arrays(:,1), arrays(:,2), 'rx')

