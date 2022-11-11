% white Gaussian signal
y = randn(1024,1);
% OR real data
%y = audioread('beepbeep.wav');
% sampling frequency in hertz
f=48000;
% microphone distance in meter
d = 0.175;
% speed of sound in meter/second
v = 343;
% maximum TDOA in samples
dn_max = round(f*d/v);
% simulated recording
x=nktp_sim(f,y,[2, 3],ones(1,8),10*ones(1,8),0);
% signal of the first array
X1 = x(:,1);
X2 = x(:,2);
% size of the data
N = size(x,1);
% calculate angle of arrival
theta = AOA_GCC(X1, X2, dn_max, N);

%% plotting
% calculate cross-correlation
rxy = genxcorr(X1, X2, N);
% shift zero lag to the middle
rxy = ifftshift(rxy);
% take middle part out, with index correction
rxy = rxy(ceil(N/2)-dn_max+1:ceil(N/2)+dn_max+1);
% time axis in relevantrange
t = -dn_max:dn_max;
% plot center part of the cross correlation
plot(t, rxy)
xlabel('$n$','Interpreter','LaTex')
ylabel('$\hat{r}_{xy}$','Interpreter','LaTex')
% matlab2tikz('noise_GCC.tikz');