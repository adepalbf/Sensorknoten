% filename: plotDOA_est_once.m 
% 
% based on: TestDOAplots.m
% ===========================================
% Aug 23, 2011 by QL
% For DOA plots, try using the [20, 20, 15]dB,  [35.11 50.15  80.31] doa, 
% [0: 0.5: 180] DOAscan
%
% For MSE-CRB: use SNR: -10 ---20 dB, source power fixed at [20 10]dB, 
% Truth DOA [35.11 50.15], [0: 0.2: 180] DOAscan, SNR defintion: mean of
% signal powers / ratio // noise_powers...
%
% =======================================
% Target: test Grid-Free searching fminsearch and compare with + RELAX
% algoss
%
% based on:  test_SAMseries_MSE_minSNR.m
%
% --------------------------------
% Updated SNR defintion 
% use the Habti's New suggestions on SNR definitions: the minimum SNR
%
%
%
%% =========== Important Parameter list ===========
% use only 2 sources for now, change locations to off-grid greatly
% use M = 12 snesors
% about: SNR = 15 dB
% about: Snapshot # = 16... etc
% The function handles independent and coherent sources...
% use 0.5 step size scan, 0: 0.5: 180
% Use one off-grid sources,  DOA = [35.21 80.31]
% % use widely separated DOA for easy DAS initialization 
 


% ============================ control running parameters ================

Flag_save_fig = 0; % 1 save figs...
Cohr_flag = 0; %0; % 1 for coherent sources, 0 for independent sources
SNR = 5; % 15; %25;  %; dB

% limited by grids
M         = 12; %10; % # of sensors
t_samples = 16; %120; %16; % 200; % # of snapshots obtained by the ULA

% % limited by measurements
% M         = 4; % very small aperture
% t_samples = 6; % very few snapshots

if ~Cohr_flag
    figpath = ['resDOAAug/DOAest_indp_M' num2str(M) '_N' num2str(t_samples) '_' num2str(SNR) 'dB/'  ];
else
    figpath = ['resDOAAug/DOAest_cohr_M' num2str(M) '_N' num2str(t_samples) '_' num2str(SNR) 'dB/'  ];
end

if ~exist(figpath, 'dir') && Flag_save_fig
    mkdir(figpath);
end

% % % const modulus noise
% % noisenew = exp(1j * 2* pi * rand(M,t_samples)); 

% Fixed Source powers
PowerDOAdB = [20; 10]; %[20 15]; %[20; 20; 15]; % in dB
PowerDOA = 10.^(PowerDOAdB/10);
amplitudeDOA = sqrt(PowerDOA);


% complex Gaussian Noise
noisePowerdB = mean(PowerDOAdB(:)) - SNR; 
noisenew = (randn(M,t_samples) + 1j* randn(M, t_samples))/sqrt(2); % noise
noisePower = 10^(noisePowerdB /10);
noisenew = noisenew * sqrt(noisePower);
% ======================================================================

Dist = ones(1, M-1); % inter-element spacing of sensors
DistTmp = cumsum([0 Dist]); % locations of the M sensors

DOAscan = 0: 0.2: 180; %0: 0.5 :180; % all possible DOA angles
% because we use cos instead of sin, so 10, 40, 55 degree -> 35, 50, 80
DOA = [35.11 50.15]; %[35.11 50.15  80.31]; % true DOA angles, on gird case 
% DOA = [95.21 120.31];
DOA = sort(DOA, 'ascend'); % must be in accend order to work right

source_no = length(DOA); % # of sources


Areal = exp(1j*pi*DistTmp' * cos(DOA*pi/180) ); % real steering vector matrix
if ~Cohr_flag % indp sources
    waveform = exp(1j*2*pi*rand(source_no, t_samples)) .* repmat(amplitudeDOA, 1, t_samples);  
else % coherent sources
    waveform = exp(1j*2*pi*rand(source_no-1, t_samples));
    waveform = [waveform;  waveform(1, :)  ];
    waveform = waveform.* repmat(amplitudeDOA , 1, t_samples);
end

y_noisefree = Areal *  waveform; % ideal noiseless measurements
y_noisy      = y_noisefree + noisenew; % noisy measurements

% steering vector matrix w.r.t all possible scanning DOA's
A = exp(1j*pi*DistTmp' * cos(DOAscan*pi/180) ); 
modulus_hat_das  = sum(abs(A'*y_noisy/M), 2 )/t_samples;


%%%%%%%%%%%%%%%%%%%%%%%%
[~, ~, p_vec_das, normal_das, ~]=fun_DASRes(y_noisy, A, modulus_hat_das,DOAscan,DOA);
if abs(normal_das) < 2*eps
    warning('DAS Abnormal!');
end
h_das = figure; 
plot(DOAscan, 10 * log10(p_vec_das + eps));
hold on;
plot(DOA,  20*log10(amplitudeDOA.'), 'ro', 'MarkerSize',10, 'LineWidth',2);
% legend('Estimates', 'Truth');
xlabel('Direction of Arrival (\circ)');
ylabel('Power (dB)');
title(['DAS ' num2str(SNR) 'dB']);
xlim([min(DOAscan) max(DOAscan)]);
ylim([-30 25]);
myboldify_frameOnly;

yminmax = get(gca, 'ylim');
line_x = repmat(DOA, 2, 1);
line_y = repmat(yminmax.', 1, length(DOA));
line(line_x, line_y, 'linewidth', 1, 'color', 'r', 'LineStyle','--'  );

if Flag_save_fig
    saveas(h_das, [figpath '/DAS_DOA_est.fig' ]);
    print('-depsc2',[figpath '/DAS_DOA_est.eps']);
end

% return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, ~, p_vec, normal]=fun_IAAEst(y_noisy, A, modulus_hat_das,DOAscan,DOA);
if abs(normal) < 2*eps
    warning('IAA Abnormal!');
end
h_iaa = figure; 
plot(DOAscan, 10 * log10(p_vec + eps));
hold on;
plot(DOA,  20*log10(amplitudeDOA.'), 'ro', 'MarkerSize',10, 'LineWidth',2);
xlabel('Direction of Arrival (\circ)');
ylabel('Power (dB)');
title(['IAA ' num2str(SNR) 'dB']);
xlim([min(DOAscan) max(DOAscan)]);
ylim([-30 25]);
myboldify_frameOnly;

yminmax = get(gca, 'ylim');
line_x = repmat(DOA, 2, 1);
line_y = repmat(yminmax.', 1, length(DOA));
line(line_x, line_y, 'linewidth', 1, 'color', 'r', 'LineStyle','--'  );

if Flag_save_fig
    saveas(h_iaa, [figpath '/IAA_DOA_est.fig' ]);
    print('-depsc2',[figpath '/IAA_DOA_est.eps']);
end



%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Detected_powers_sam0, Distance_sam0, p_vec_sam0, normal_sam0, NoisePower_sam0] = ...
    fun_SAM0Res(y_noisy, A, modulus_hat_das,DOAscan,DOA);
if abs(normal_sam0) < 2*eps
    warning('SAM0 Abnormal!');
end
h_sam0 = figure; 
plot(DOAscan, 10 * log10(p_vec_sam0 + eps));
hold on;
plot(DOA,  20*log10(amplitudeDOA.'), 'ro');
legend('Estimates', 'Truth');
xlabel('Direction of Arrival (\circ)');
ylabel('Power (dB)');
title(['SAMV-0 ' num2str(SNR) 'dB']);
xlim([min(DOAscan) max(DOAscan)]);
myboldify_frameOnly;
if Flag_save_fig
   saveas(h_sam0, [figpath '/SAM0_DOA_est.fig' ]);
   print('-depsc2',[figpath '/SAM0_DOA_est.eps']);
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Detected_powers_sam1, Distance_sam1, p_vec_sam1, normal_sam1, NoisePower_sam1] = ...
    fun_SAM1Res(y_noisy, A, modulus_hat_das,DOAscan,DOA);
if abs(normal_sam1) < 2*eps
    warning('SAM1 Abnormal!');
end
h_sam1 = figure; 
plot(DOAscan, 10 * log10(p_vec_sam1 + eps));
hold on;
plot(DOA,  20*log10(amplitudeDOA.'), 'ro', 'MarkerSize',10, 'LineWidth',2);
% legend('Estimates', 'Truth');
xlabel('Direction of Arrival (\circ)');
ylabel('Power (dB)');
title(['SAMV-1 ' num2str(SNR) 'dB']);
xlim([min(DOAscan) max(DOAscan)]);
ylim([-30 25]);
myboldify_frameOnly;

yminmax = get(gca, 'ylim');
line_x = repmat(DOA, 2, 1);
line_y = repmat(yminmax.', 1, length(DOA));
line(line_x, line_y, 'linewidth', 1, 'color', 'r', 'LineStyle','--'  );

if Flag_save_fig
   saveas(h_sam1, [figpath '/SAM1_DOA_est.fig' ]);
   print('-depsc2',[figpath '/SAM1_DOA_est.eps']);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Detected_powers_sam3, Distance_sam3, p_vec_sam3, normal_sam3, NoisePower_sam3] = ...
    fun_SAM3Res(y_noisy, A, modulus_hat_das,DOAscan,DOA);
if abs(normal_sam3) < 2*eps
    warning('SAM3 Abnormal!');
end
h_sam3 = figure; 
plot(DOAscan, 10 * log10(p_vec_sam3 + eps));
hold on;
plot(DOA,  20*log10(amplitudeDOA.'), 'ro', 'MarkerSize',10, 'LineWidth',2);
% legend('Estimates', 'Truth');
xlabel('Direction of Arrival (\circ)');
ylabel('Power (dB)');
title(['SAMV-3 ' num2str(SNR) 'dB']);
xlim([min(DOAscan) max(DOAscan)]);
ylim([-30 25]);
myboldify_frameOnly;

yminmax = get(gca, 'ylim');
line_x = repmat(DOA, 2, 1);
line_y = repmat(yminmax.', 1, length(DOA));
line(line_x, line_y, 'linewidth', 1, 'color', 'r', 'LineStyle','--'  );

if Flag_save_fig
   saveas(h_sam3, [figpath '/SAM3_DOA_est.fig' ]);
   print('-depsc2',[figpath '/SAM3_DOA_est.eps']);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Detected_powers, Distance, p_vec, normal, noisepower] = fun_SPICEplusRes(y_noisy, A, modulus_hat_das,DOAscan,DOA);
if abs(normal) < 2*eps
    warning('SPICEplus Abnormal!');
end
h_spice = figure; 
plot(DOAscan, 10 * log10(p_vec + eps));
hold on;
plot(DOA, 20*log10(amplitudeDOA.'), 'ro', 'MarkerSize',10, 'LineWidth',2);
xlabel('Direction of Arrival (\circ)');
ylabel('Power (dB)');
title(['SPICE+ ' num2str(SNR) 'dB']);
xlim([min(DOAscan) max(DOAscan)]);
ylim([-30 25]);
myboldify_frameOnly;

yminmax = get(gca, 'ylim');
line_x = repmat(DOA, 2, 1);
line_y = repmat(yminmax.', 1, length(DOA));
line(line_x, line_y, 'linewidth', 1, 'color', 'r', 'LineStyle','--'  );

if Flag_save_fig
    saveas(h_spice, [figpath '/SPICEplus_DOA_est.fig' ]);
    print('-depsc2',[figpath '/SPICEplus_DOA_est.eps']);
end







%% ===   Res1_ML implementation.. 
tic;
% Must compute the initialziation values here

% [Detected_powers_sam3, Distance_sam3, p_vec_sam3, normal_sam3, NoisePower_sam3] = ...
%     fun_SAM3Res(y_noisy, A, modulus_hat_das,DOAscan,DOA);
normal_in = normal_sam3;
initDOA = DOA + Distance_sam3;
initNoisePower = NoisePower_sam3;
initPower = Detected_powers_sam3;


[p_k_mat_r1ML,  noisePower_r1ML, theta_k_mat_r1ML, Distance_r1ML, normal_r1ML]=...
    fun_Res1_MLRes(y_noisy, A, DOA, initPower, initDOA, initNoisePower, normal_in);
time = toc;
disp(['Res1-ML algo running time ==' num2str(time) ' sec.' ]);

if abs(normal_r1ML) < 2*eps
    warning('Res1ML Abnormal!');
end

h_r1ML = figure; 
plot(theta_k_mat_r1ML, 10 * log10(p_k_mat_r1ML + eps), 'b*', 'MarkerSize',10, 'LineWidth',2);
hold on;
plot(DOA,  20*log10(amplitudeDOA.'), 'ro', 'MarkerSize',10, 'LineWidth',2);
% legend('Estimates', 'Truth');
xlabel('Direction of Arrival (\circ)');
ylabel('Power (dB)');
title(['Res1 ML ' num2str(SNR) 'dB']);
xlim([min(DOAscan) max(DOAscan)]);
ylim([-30 25]);
myboldify_frameOnly;

yminmax = get(gca, 'ylim');
line_x = repmat(DOA, 2, 1);
line_y = repmat(yminmax.', 1, length(DOA));
line(line_x, line_y, 'linewidth', 1, 'color', 'r', 'LineStyle','--'  );

if Flag_save_fig
   saveas(h_r1ML, [figpath '/r1ML_DOA_est.fig' ]);
   print('-depsc2',[figpath '/r1ML_DOA_est.eps']);
end





%% ===   SAMV_k_ML implementation.. 
Distance_samkML = {};
for samk = [  1 3]
    tic;
    % Must compute the initialziation values here
    normal_in = normal_sam3;
    initDOA = DOA + Distance_sam3;
    initNoisePower = NoisePower_sam3;
    initPower = Detected_powers_sam3;


    [p_k_mat,  noisePower, theta_k_mat, Distance, normal]=...
        fun_SAMk_MLRes(y_noisy, A, DOA, initPower, initDOA, initNoisePower, normal_in, samk);
    time = toc;
    disp(['SAMV-'  num2str(samk) '-ML algo running time ==' num2str(time) ' sec.' ]);
    Distance_samkML{end+1} = Distance;

    if abs(normal ) < 2*eps
        warning(['SAMV-' num2str(samk) '-ML Abnormal!']);
    end

    h_ML = figure; 
    plot(theta_k_mat, 10 * log10(p_k_mat + eps), 'b*', 'MarkerSize',10, 'LineWidth',2);
    hold on;
    plot(DOA,  20*log10(amplitudeDOA.'), 'ro', 'MarkerSize',10, 'LineWidth',2);
%     legend('Estimates', 'Truth');
    xlabel('Direction of Arrival (\circ)');
    ylabel('Power (dB)');
    title(['SAMV-' num2str(samk) '-ML ' num2str(SNR) 'dB']);
    xlim([min(DOAscan) max(DOAscan)]);
    ylim([-30 25]);
    myboldify_frameOnly;
    
    yminmax = get(gca, 'ylim');
    line_x = repmat(DOA, 2, 1);
    line_y = repmat(yminmax.', 1, length(DOA));
    line(line_x, line_y, 'linewidth', 1, 'color', 'r', 'LineStyle','--'  );

    if Flag_save_fig
       saveas(h_ML, [figpath '/SAMV-' num2str(samk) '-ML_DOA_est.fig' ]);
       print('-depsc2',[figpath '/SAMV-' num2str(samk) '-ML_DOA_est.eps']);
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Distance, p_vec, normal] = ...
[pseu_p_select_music, theta_k_select_music, p_vec, normal, ~] = ...    
    fun_MUSIC_GridFree_Res(y_noisy, A, modulus_hat_das,DOAscan,DOA);
if abs(normal) < 2*eps
    warning('MUSIC Abnormal!');
end
h_music = figure; 
plot(DOAscan, 10 * log10(p_vec + eps));
hold on;
xlim([min(DOAscan) max(DOAscan)]);
xlabel('Direction of Arrival (\circ)');
ylabel('Pseudo Power / dB');
title(['MUSIC ' num2str(SNR) 'dB']);
myboldify_frameOnly;


%     %--- this is old-fashioned, the amplitude is not relevant...
%     plot(DOA,  20*log10([10 10 3]), 'ro');
% ----- Using lines instead of the red circles ...
yminmax = get(gca, 'ylim');
line_x = repmat(DOA, 2, 1);
line_y = repmat(yminmax.', 1, length(DOA));
line(line_x, line_y, 'linewidth', 1, 'color', 'r', 'LineStyle','--'  );
if Flag_save_fig
    saveas(h_music, [figpath '/MUSIC_DOA_est.fig' ]);
    print('-depsc2',[figpath '/MUSIC_DOA_est.eps']);
end






