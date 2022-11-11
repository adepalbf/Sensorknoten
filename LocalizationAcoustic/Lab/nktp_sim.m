function x=nktp_sim(FREQ, SPEECH, LOC, G, SNR, NO_NOISE, arrays, d )

% x=nktp_sim(FREQ, SPEECH, LOC, G, SNR, NO_NOISE )
%  generate simulation data.
% inputs:
%  1- FREQ: the sampling frequency.
%  2- SPEECH: the speech signal.
%  3- LOC: the true location of the acoustical source.
%  4- G: the gain at the 8 microphones.
%  5- SNR: the SNR at the 8 microphones.
%  6- NO_NOISE: 1 for non noisy signals and 0 for noise signals.
% output
%  x: the output signal at the microphones.


% if nargin ~= 6, error('Please provide the 6 parameters!'), end

% experiment setup
% arrays = [2.69, 0.19; 
%            5.19, 3.21;
%            2.70, 5.27;
%            0.37, 3.54];
% d  = 0.175;
vs = 343;
mics(1, :) = [arrays(1,1)  - d/2, arrays(1,2)      ];
mics(2, :) = [arrays(1,1)  + d/2, arrays(1,2)      ];
mics(3, :) = [arrays(2,1)       , arrays(2,2) - d/2];
mics(4, :) = [arrays(2,1)       , arrays(2,2) + d/2];
mics(5, :) = [arrays(3,1)  + d/2, arrays(3,2)      ];
mics(6, :) = [arrays(3,1)  - d/2, arrays(3,2)      ];
mics(7, :) = [arrays(4,1)       , arrays(4,2) + d/2];
mics(8, :) = [arrays(4,1)       , arrays(4,2) - d/2];

% calculate delays
dist2mics = zeros(8,1);
for k=1:8
    p2   = mics(k, :);
    xy_dif = LOC - p2;
    dist2mics(k) = sqrt(xy_dif*xy_dif');
end
dist_dif = dist2mics(2:2:end) - dist2mics(1:2:end-1);
delay_sec = dist_dif / vs;
n_delay = round(delay_sec * FREQ);

% generate non noisy signals
max_n_delay = max(abs(n_delay));
N = length(SPEECH) - max_n_delay;
x = zeros(N, 8);
for k=1:4
    n = n_delay(k);
    if n > 0
        x(:, 2*k-1) = G(2*k-1) * SPEECH(n+1:n+N);
        x(:, 2*k)   = G(2*k)   * SPEECH(1:N);
    else
        x(:, 2*k-1) = G(2*k-1) * SPEECH(1:N);
        x(:, 2*k)   = G(2*k)   * SPEECH(-n+1:-n+N);
    end    
end

% add noise
if ~NO_NOISE
    for i=1:8
        sig_p = sum(abs(x(:, i)).^2)/length(x(:, i));
        var_noise = sig_p/(10^(SNR(i)/20));
        noise = randn(size(x(:, i)));
        noise_p = sum(abs(noise).^2)/length(noise);
        noise = var_noise*noise/noise_p;
        x(:, i) = x(:, i) + noise;
    end
end

end
