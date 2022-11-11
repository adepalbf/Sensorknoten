% ============================ control running parameters ================
MC = 1; % # of plotting together...
Flag_save_fig = 0; % 1 save figs...
Cohr_flag = 0; %0; % 1 for coherent sources, 0 for independent sources

M = 2;

% speed of sound in meter/second
v = 343;

%% ======================================================================
% Array
dx = 0.15;
dy = 0.15;
Nx = 2;
Ny = 2;

b_no_noise = 0;
SNR = 10; % 15; %25;  %; dB
arrays = [5,0,0;10,5,0;5,10,0;0,5,0];
direction = [0;pi/2;pi;3*pi/2];
N = length(direction);
Arrays = array.empty(N,0);

h_arr = figure;
for ii = 1:N
    Arrays(ii) = array(arrays(ii,:), direction(ii), dx, dy, Nx, Ny);
    Arrays(ii).plot_array()
    text(Arrays(ii).pos(1)+0.3,Arrays(ii).pos(2),"Node "+ii,'VerticalAlignment','bottom','HorizontalAlignment','left')

    hold on
end

%%
phi_scan = (-90:0.5:90)*pi/180; %0: 0.2: 180; % all possible DOA angles

%% ======================================================================

phi_source = 50*pi/180; % Azimuth/Azimut
theta_source = 90*pi/180; % Inclination/Polwinkel
r_source = 4;
% phi_source = [30*pi/180, 45*pi/180];
% theta_source = [60*pi/180, 90*pi/180];
% r_source = [1, 2];
pos_source = r_source.*[sin(theta_source).*cos(phi_source); sin(theta_source).*sin(phi_source); cos(theta_source)];

% Fixed Source powers
PowerDOAdB = 25; %[5; 3; 4]; %[20; 20; 15]; % in dB
PowerDOA = 10.^(PowerDOAdB/10);
amplitudeDOA = sqrt(PowerDOA);

%% ======================================================================
var_data = 3; % 1 - sine, 2 - noise, 3 - beepbeep1
T = 10;
f = 1000;
fs = 48000;
t = (1:T*fs)/fs;
t_samples = T*fs; 
switch var_data
    case 1
        % sine
        y = sin(2*pi*f*t);
    case 2
        % white Gaussian signal
        y = randn(T*fs,1);
    case 3
        % beepbeep
        [y, Fs] = audioread('beepbeep.wav');  
        y = resample(y,fs,Fs);
        y = repmat(y, ceil(t_samples/length(y)));
        y = y(1:t_samples);
end


%% ======================================================================
if ~Cohr_flag
    figpath = ['resDOASep4/DOAest_indp_M' num2str(M) '_N' num2str(t_samples) '_' num2str(SNR) 'dB/'  ];
else
    figpath = ['resDOASep4/DOAest_cohr_M' num2str(M) '_N' num2str(t_samples) '_' num2str(SNR) 'dB/'  ];
end

if ~exist(figpath, 'dir') && Flag_save_fig
    mkdir(figpath);
end
%%
plot(pos_source(1),pos_source(2),'o', 'MarkerSize', 10)
text(pos_source(1),pos_source(2),'Real position','VerticalAlignment','top','HorizontalAlignment','right')

hold on
%%
beta = zeros(N,1);
P = zeros(N,1);
beta_pos = zeros(N,3);

for ind = 1:MC
    
    
    for ii = 1:N
        % simulated recording
        y_MC = Arrays(ii).record( fs,y,pos_source,amplitudeDOA, SNR, b_no_noise);
        
        % estimate AOA
        [beta(ii), P(ii)] = Arrays(ii).estimate_AOA(y_MC, fs,'AOA_CC');
        
        figure(h_arr)
        % debug code
         beta_pos(ii,:) = Arrays(ii).global_coords(beta(ii), r_source);
         error = Arrays(ii).error(beta(ii), phi_source, r_source);
         fprintf('Error(%d) = %.2f°\n', ii,error*180/pi)
    end
    %%
    
    % transform theta to global coordinates
    theta = (beta + direction);
    % calculate least-squares estimate of the source position
    A = [tan(theta)./sqrt(tan(theta).^2+1), -1./sqrt(tan(theta).^2+1)];
    b = arrays(:,2) - tan(theta).*arrays(:,1);
    b = b./sqrt(tan(theta).^2+1);
    p_hat = -(A'*A)\A'*b;
    
    plot(p_hat(1),p_hat(2), 'x', 'MarkerSize', 10);
    text(p_hat(1),p_hat(2),'Estimated position','VerticalAlignment','bottom','HorizontalAlignment','left')
    fprintf('Error_abs = %.1fm\n', sqrt(sum(p_hat-pos_source(1:2)).^2))

end
%%
% figure
% plot(phi_source*180/pi*[1,1],  [0, PowerDOAdB], 'r', 'MarkerSize',10, 'LineWidth',2);
% % legend('Estimates', 'Truth');
% xlabel('Direction of Arrival (\circ)');
% ylabel('Power (dB)');
% if Cohr_flag
%     title(['Coherent GCC SNR ' num2str(SNR) 'dB']);
% else
%     title(['Independent GCC SNR ' num2str(SNR) 'dB']);
% end
% xlim([min(phi_scan) max(phi_scan)]*180/pi);
% ylim([0 PowerDOAdB]);
% 
% % myboldify_frameOnly;
% if Flag_save_fig
%     saveas(h_gcc, [figpath '/GCC_DOA_est.fig' ]);
%     print('-depsc2',[figpath '/GCC_DOA_est.eps']);
% end



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
