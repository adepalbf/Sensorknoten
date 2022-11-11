function x=nktp_sim(fs, y, pos_source, gain_source, SNR, b_no_noise, mics)

% x=nktp_sim(FREQ, SPEECH, LOC, G, SNR, NO_NOISE )
%  generate simulation data.
% inputs:
%  1- FREQ: the sampling frequency.
%  2- SPEECH: the speech signal.
%  3- LOC: the true location of the acoustical source.
%  4- G: the gain at each microphone.
%  5- SNR: the SNR at each microphone.
%  6- NO_NOISE: 1 for non noisy signals and 0 for noise signals.
% output
%  x: the output signal at the microphones.

vs = 343;

M = size(mics,2); % # of sensors
N_source = size(pos_source, 2); % # of sources
% calculate delays
dist2mics = zeros(N_source, M);
for kk=1:M
    for ll=1:N_source
        xy_dif = pos_source(:,ll) - mics(:, kk);
        dist2mics(ll, kk) = sqrt(xy_dif'*xy_dif);
    end
end

tau = dist2mics / vs;
tau = tau - min(tau);
% legendStr = {};
x = zeros(length(y), M);
for kk=1:M
    sys = thiran(tau(kk), 1/fs);
    x(:,kk) = lsim(sys,y);
    
%     x_temp = lsim(sys,y);
%     x(:,kk) = gain_source(kk) * x_temp;

%     legendStr{end+1} = sprintf('Mic %d delay %d', kk, tau(kk));
end

x = gain_source .* x;

% n_delay = round(tau * fs);
% % generate non noisy signals
% max_n_delay = max(abs(n_delay), [], 'all');
% N = length(y);% - max_n_delay;
% y = [zeros(max_n_delay+1, 1); y];
%
% x2 = zeros(N, M);
%
% for kk=1:M
%     n = max_n_delay-n_delay(kk);
%     x2(1:N, kk) = gain_source(kk) * y(n+(1:N));
%     legendStr{end+1} = sprintf('Mic %d delay %d', kk, n_delay(kk));
% end
% figure
% plot(x)
% hold all
% ax = gca;
% ax.ColorOrderIndex = 1;
% plot(x2, '--')
% legend(legendStr)

% add noise
if ~b_no_noise
        Px = sum(x.^2)/size(x,1);
        var_noise = Px/(10^(SNR/20));
        noise = var_noise.*randn(size(x));
%         Pn = sum(noise.^2)/size(noise,1);
%         SNR_eff = 10*log10(Px/Pn);
        x = x + noise;
end

end
