% ============================ control running parameters ================
Flag_save_fig = 0; % 1 save figs...
Cohr_flag = 0   ; %0; % 1 for coherent sources, 0 for independent sources
SNR = 25; % 15; %25;  %; dB

% white Gaussian signal
y = randn(2^14,1);
% OR real data
% y = audioread('beepbeep.wav');
y = y / rms(y);
% sampling frequency in hertz
f=48000;
% microphone distance in meter
d = 0.175;
% speed of sound in meter/second
v = 343;
% maximum TDOA in samples
dn_max = round(f*d/v);

b_noisy = 1;

% M = 8;
% Dist = ones(1, M-1); % inter-element spacing of sensors
% DistX = cumsum([0 Dist]); % locations of the M sensors
% DistY = zeros(size(DistX));
% arrays = [DistX',DistY']*d; % sensors in rows

arrays = [2.69, 0.19;
    5.19, 3.21;
    2.70, 5.27;
    0.37, 3.54];

mics(1, :) = [arrays(1,1)  - d/2, arrays(1,2)      ];
mics(2, :) = [arrays(1,1)  + d/2, arrays(1,2)      ];
mics(3, :) = [arrays(2,1)       , arrays(2,2) - d/2];
mics(4, :) = [arrays(2,1)       , arrays(2,2) + d/2];
mics(5, :) = [arrays(3,1)  + d/2, arrays(3,2)      ];
mics(6, :) = [arrays(3,1)  - d/2, arrays(3,2)      ];
mics(7, :) = [arrays(4,1)       , arrays(4,2) + d/2];
mics(8, :) = [arrays(4,1)       , arrays(4,2) - d/2];

M = size(arrays,1); % # of sensors

DOAscan = -90: 0.5 :90; %0: 0.2: 180; % all possible DOA angles

% limited by grids
t_samples = 120; %16; % 200; % # of snapshots obtained by the ULA

if ~Cohr_flag
    figpath = ['resDOASep4/DOAest_indp_M' num2str(M) '_N' num2str(t_samples) '_' num2str(SNR) 'dB/'  ];
else
    figpath = ['resDOASep4/DOAest_cohr_M' num2str(M) '_N' num2str(t_samples) '_' num2str(SNR) 'dB/'  ];
end

if ~exist(figpath, 'dir') && Flag_save_fig
    mkdir(figpath);
end

% % % const modulus noise
% % noisenew = exp(1j * 2* pi * rand(M,t_samples));
pos = [8.19, 5.74];
DOA = atan(pos(1)/pos(2))*180/pi;
DOA = sort(DOA, 'ascend'); % must be in accend order to work right
source_no = length(DOA); % # of sources

% Fixed Source powers
PowerDOAdB = 25; %[5; 3; 4]; %[20; 20; 15]; % in dB
PowerDOA = 10.^(PowerDOAdB/10);
amplitudeDOA = sqrt(PowerDOA);

% ======================================================================

k = 1*pi/d*[cos(DOA*pi/180); sin(DOA*pi/180)]; % sources in columns

% Areal = exp(1i*pi*DistTmp' * cos(DOA*pi/180) ); % real steering vector matrix
Areal = exp(1i*mics*k);

% steering vector matrix w.r.t all possible scanning DOA's
k = 1*pi/d*[cos(DOAscan*pi/180); sin(DOAscan*pi/180)]; % sources in columns
A = exp(1i*mics*k);
% A = exp(1i*pi*DistTmp' * cos(DOAscan*pi/180) );
%%
% simulated recording
y_noisy = nktp_sim(f,y/rms(y),pos,amplitudeDOA*ones(1,2*M),SNR*ones(1,2*M),b_noisy, mics)';

modulus_hat_das  = sum(abs(A'*y_noisy/M), 2 )/t_samples;


%%
h_sam0 = figure;
[p_vec_sam0, NoisePower_sam0] = SAMV0(y_noisy, A, modulus_hat_das);
[Detected_powers_sam0, Distance_sam0, normal_sam0]=SAMV0_result(p_vec_sam0, DOAscan, DOA);

if abs(normal_sam0) < 2*eps
    warning('SAM0 Abnormal!');
end
figure(h_sam0)
plot(DOAscan, 10 * log10(p_vec_sam0 + eps));
hold on;


plot(DOA,  PowerDOAdB, 'ro', 'MarkerSize',10, 'LineWidth',2);
% legend('Estimates', 'Truth');
xlabel('Direction of Arrival (\circ)');
ylabel('Power (dB)');
if Cohr_flag
    title(['Coherent SAMV-0 ' num2str(SNR) 'dB']);
else
    title(['Independent SAMV-0 ' num2str(SNR) 'dB']);
end
% xlim([min(DOAscan) max(DOAscan)]);
% ylim([-35 10]);

% myboldify_frameOnly;
if Flag_save_fig
    saveas(h_sam0, [figpath '/SAM0_DOA_est.fig' ]);
    print('-depsc2',[figpath '/SAM0_DOA_est.eps']);
end

%%

h_gcc = figure;
alpha = 0.8;
% simulated recording
I = zeros(t_samples,M-1);
N = 1024;
[p_hat, theta, rxy] = estimate_position_AOA_GCC(y_noisy, dn_max, N, arrays);

theta = theta * 180/pi;
[Distance_gcc, normal_gcc]=GCC_result(theta, DOA);

theta = mean(theta);
Detected_powers = mean(Detected_powers);

if abs(normal_gcc) < 2*eps
    warning('GCC Abnormal!');
end

plot(theta,  20*log10(Detected_powers.'), 'kx', 'MarkerSize',10, 'LineWidth',2);
hold on;

x = cell(size(rxy,2), 1);
rxy_i = cell(size(rxy,2), 1);
rxy_i_dB = cell(size(rxy,2), 1);
for ii = 1:size(rxy,2)
    rxy_shift = ifftshift(rxy(:,ii));
    range = -dn_max(ii):dn_max(ii);
    rxy_i{ii} = rxy_shift(t_samples/2+range);
    
    rxy_i_dB{ii} = 20 * log10(abs(rxy_i{ii}) + eps);
    x{ii} = asin(range/dn_max(ii))*180/pi;
    plot(x{ii}, rxy_i_dB{ii});
    hold on;
end
plot(DOA,  PowerDOAdB, 'ro', 'MarkerSize',10, 'LineWidth',2);
% legend('Estimates', 'Truth');
xlabel('Direction of Arrival (\circ)');
ylabel('Power (dB)');
if Cohr_flag
    title(['Coherent GCC ' num2str(SNR) 'dB']);
else
    title(['Independent GCC ' num2str(SNR) 'dB']);
end
% xlim([min(DOAscan) max(DOAscan)]);
% ylim([-30 0]);

% myboldify_frameOnly;
if Flag_save_fig
    saveas(h_gcc, [figpath '/GCC_DOA_est.fig' ]);
    print('-depsc2',[figpath '/GCC_DOA_est.eps']);
end
