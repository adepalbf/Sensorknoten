% ============================ control running parameters ================
MC = 10; % # of plotting together...
Flag_save_fig = 0; % 1 save figs...
Cohr_flag = 0; %0; % 1 for coherent sources, 0 for independent sources
SNR = 10; % 15; %25;  %; dB
t_samples = 1000; 
fs = 48000;

M = 2;
b_no_noise = 1;

if ~Cohr_flag
    figpath = ['resDOASep4/DOAest_indp_M' num2str(M) '_N' num2str(t_samples) '_' num2str(SNR) 'dB/'  ];
else
    figpath = ['resDOASep4/DOAest_cohr_M' num2str(M) '_N' num2str(t_samples) '_' num2str(SNR) 'dB/'  ];
end

if ~exist(figpath, 'dir') && Flag_save_fig
    mkdir(figpath);
end


% speed of sound in meter/second
v = 343;


%% ======================================================================
% Array

D_max = 0.19;

% microphone distance in meter
d = D_max/(M-1);%0.175;
u0 = [0;1;0];
ULA = (0:M-1)*d-D_max/2;
ULA = u0*ULA;
phi_T = 0;
theta_T = 0;
T = [cosd(phi_T), -sind(phi_T), 0; sind(phi_T), cosd(phi_T) 0; 0 0 1];
arrays = T*ULA; % sensors in rows

phi_scan = (-90:0.5:90)*pi/180; %0: 0.2: 180; % all possible DOA angles


% maximum TDOA in samples
difference = arrays(:,2:end)-arrays(:,1);
distance = diag(sqrt(difference'*difference));
dn_max = round(fs*distance/v);

phi_resolve = asin(1./dn_max)*180/pi;

%% ======================================================================

phi_source = 10; % Azimuth/Azimut
theta_source = 90; % Inclination/Polwinkel
r_source = 10;
% phi_source = [30; 45]*pi/180;
% r_source = [1; 2];
pos_source = r_source.*[sind(theta_source)*cosd(phi_source); sind(theta_source)*sind(phi_source); cosd(theta_source)];

% Fixed Source powers
PowerDOAdB = 25; %[5; 3; 4]; %[20; 20; 15]; % in dB
PowerDOA = 10.^(PowerDOAdB/10);
amplitudeDOA = sqrt(PowerDOA);

%% ======================================================================
var_data = 3; % 1 - sine, 2 - noise, 3 - beepbeep
T = 10;
f = 1000;
y = audio_data(var_data, fs, T, f);

%% ======================================================================

h_gcc = figure;

for ind = 1:MC
    % simulated recording
    y_MC = nktp_sim(fs,y,pos_source,amplitudeDOA*ones(1,M),SNR*ones(1,M),b_no_noise, arrays);
    [beta, Detected_powers] = AOA_CC(y_MC, dn_max);

%     rxy_av = rxy_av+rxy/MC;
    beta = beta * 180/pi;
%     [Distance_gcc, normal_gcc]=GCC_result(beta, phi_source*180/pi);
    
    theta = mean(beta);
    P = mean(Detected_powers);
    
%     if abs(normal_gcc) < 2*eps
%         warning('GCC Abnormal!');
%     end
    figure(h_gcc)
    stem(beta,  10*log10(Detected_powers.'));
    hold all;
    stem(theta,  10*log10(P), 'k', 'MarkerSize',10, 'LineWidth',2);
end

% x = cell(size(rxy,2), 1);
% rxy_i = cell(size(rxy,2), 1);
% rxy_i_dB = cell(size(rxy,2), 1);
% for ii = 1:size(rxy,2)
%     rxy_shift = ifftshift(rxy_av(:,ii));
%     range = -dn_max(ii):dn_max(ii);
%     rxy_i{ii} = rxy_shift(t_samples/2+range);
    
%     rxy_i_dB{ii} = 20 * log10(abs(rxy_i{ii}) + eps);
%     x{ii} = asin(range/dn_max(ii))*180/pi;
%     plot(x{ii}, rxy_i_dB{ii});
%     hold on;
% end
plot(phi_source*[1,1],  [0, PowerDOAdB], 'r', 'MarkerSize',10, 'LineWidth',2);
% legend('Estimates', 'Truth');
xlabel('Direction of Arrival (\circ)');
ylabel('Power (dB)');
if Cohr_flag
    title(['Coherent GCC SNR ' num2str(SNR) 'dB']);
else
    title(['Independent GCC SNR ' num2str(SNR) 'dB']);
end
xlim([min(phi_scan) max(phi_scan)]*180/pi);
ylim([0 PowerDOAdB]);

% myboldify_frameOnly;
if Flag_save_fig
    saveas(h_gcc, [figpath '/GCC_DOA_est.fig' ]);
    print('-depsc2',[figpath '/GCC_DOA_est.eps']);
end



%%
% [Areal, A_scan] = scanning_arrays(f, phi_source, arrays, phi_scan);

% % complex Gaussian Noise
% noisePowerdB = mean(PowerDOAdB(:)) - SNR;
% noisePower = 10^(noisePowerdB /10);
% 
% y_noisefree = cell(MC, 1);
% y_noisy = cell(MC, 1);
% modulus_hat_das = cell(MC, 1);
% for ind = 1:MC
%     
%     noisenew = (randn(M,t_samples) + 1j* randn(M, t_samples))/sqrt(2) * sqrt(noisePower); % noise
%     
%     if ~Cohr_flag % indp sources
%         waveform = exp(1j*2*pi*rand(source_no, t_samples)) .* repmat(amplitudeDOA, 1, t_samples);
%     else % coherent sources
%         waveform = exp(1j*2*pi*rand(source_no-1, t_samples));
%         waveform = [waveform;  waveform(1, :)  ]; %#ok<AGROW>
%         waveform = waveform .* repmat(amplitudeDOA , 1, t_samples);
%     end
%     
%     y_noisefree{ind} = Areal *  waveform; % ideal noiseless measurements
%     y_noisy{ind}     = y_noisefree{ind} + noisenew; % noisy measurements
%     
%     %     % simulated recording
%     %     y_noisy{ind} = nktp_sim(f,y/rms(y),pos,amplitudeDOA*ones(1,M),SNR*ones(1,M),b_noisy, arrays, t_samples)';
%     
%     modulus_hat_das{ind}  = sum(abs(A_scan'*y_noisy{ind}/M), 2 )/t_samples;
% end
% 
% %%
% h_sam0 = figure;
% for ind = 1:MC
%     [p_vec_sam0, NoisePower_sam0] = SAMV0(y_noisy{ind}, A_scan, modulus_hat_das{ind});
%     [Detected_powers_sam0, Distance_sam0, normal_sam0]=SAMV0_result(p_vec_sam0, phi_scan, phi_source);
%     
%     if abs(normal_sam0) < 2*eps
%         warning('SAM0 Abnormal!');
%     end
%     figure(h_sam0)
%     plot(phi_scan, 10 * log10(p_vec_sam0 + eps));
%     hold on;
% end
% 
% plot(phi_source,  PowerDOAdB, 'ro', 'MarkerSize',10, 'LineWidth',2);
% % legend('Estimates', 'Truth');
% xlabel('Direction of Arrival (\circ)');
% ylabel('Power (dB)');
% if Cohr_flag
%     title(['Coherent SAMV-0 ' num2str(SNR) 'dB']);
% else
%     title(['Independent SAMV-0 ' num2str(SNR) 'dB']);
% end
% % xlim([min(DOAscan) max(DOAscan)]);
% % ylim([-35 10]);

% myboldify_frameOnly;
if Flag_save_fig
    saveas(h_sam0, [figpath '/SAM0_DOA_est.fig' ]);
    print('-depsc2',[figpath '/SAM0_DOA_est.eps']);
end
