function [SqrErr, PowerSE] = ... % -10 ~ 20; % dB
    Compute_algos_StdErr(SNR_in, t_samples_in, M_in, cohr_flag) % about OK...
% updated Aug. 24, 2011 by QL
% =====  New True signal powers 
%
% new SNR definition...
% new for SAMV-ML series...
% Outputs: 
% SqrErr in DOA, in degree^2, not in dB
% PowerSE in power estimates, not in dB
%
% --------------------------------
% use the Habti mean SNR definitions.....
%
%% =========== modifications list ===========
%  grid size 0.2 deg, [0: 0.2 : 180]


algo_list = {'PER', 'IAA', 'SAMV1', 'SAMV2', 'SPICE+', 'AMV-SML', 'SAMV1-SML', 'SAMV2-SML', 'MUSIC'};
% algo_list = {'PER', 'IAA'};



% % for debugging only ....
% disp('========== Debugging mode at Function: Compute_algos_StdErr!!! ');
% SNR_in = 0; % -10 ... 25
% t_samples_in = 16; % 16 or 120
% M_in = 12;
% cohr_flag = 0;  % or 1 for coherent sources.
% % end for debugging only...





% ============================ control running parameters ================
SqrErr = cell(9, 1);
PowerSE = cell(9, 1);

M         = M_in; % # of sensors
t_samples = t_samples_in; %20;  % # of snapshots obtained by the ULA

% Fixed Source powers
PowerDOAdB = [5; 3]; % in dB
PowerDOA = 10.^(PowerDOAdB/10);
amplitudeDOA = sqrt(PowerDOA);




% complex Gaussian Noise
SNR = SNR_in; % input parameters
noisePowerdB = mean(PowerDOAdB(:)) - SNR; 
noisenew = (randn(M,t_samples) + 1j* randn(M, t_samples))/sqrt(2); % noise
noisePower = 10^(noisePowerdB /10);
noisenew = noisenew * sqrt(noisePower);




% ======================================================================


Dist = ones(1, M-1); % inter-element spacing of sensors
DistTmp = cumsum([0 Dist]); % locations of the M sensors

DOAscan = 0: 0.2 :180; % all possible DOA angles
% because we use cos instead of sin, so 10, 40, 55 degree -> 35, 50, 80
DOA = [35.11 50.15]; % true DOA angles, off gird case 
% DOA = [35 50.01 80]; % true DOA angles, on gird case 
DOA = sort(DOA, 'ascend'); % must be in accend order to work right

%% SNR       = [25 25 25].'; % must be a col vec
source_no = length(DOA); % # of sources


Areal = exp(1j*pi*DistTmp' * cos(DOA*pi/180) ); % real steering vector matrix

if ~cohr_flag % indp sources
    waveform = exp(1j*2*pi*rand(source_no, t_samples)) .* repmat(amplitudeDOA, 1, t_samples); 
else % cohr sources
    waveform = exp(1j*2*pi*rand(source_no-1, t_samples));
    waveform = [waveform;  waveform(1, :)  ];
    waveform = waveform.* repmat(amplitudeDOA , 1, t_samples);
end

    
    
y_noisefree = Areal *  waveform; % ideal noiseless measurements
y_noisy      = y_noisefree + noisenew; % noisy measurements

% steering vector matrix w.r.t all possible scanning DOA's
A = exp(1j*pi*DistTmp' * cos(DOAscan*pi/180) ); 
modulus_hat_das  = sum(abs(A'*y_noisy/M), 2 )/t_samples;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ===========  Begin Calling Functions... ==================== 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Detected_powers, Distance, ~, normal ~]=fun_DASRes(y_noisy, A, modulus_hat_das,DOAscan,DOA);
if ~normal
%     warning('DAS Abnormal!');
    SqrErr{1} = NaN;
    PowerSE{1} = NaN;
else
    SqrErr{1} = Distance * Distance';
    power_dif = Detected_powers - PowerDOA;
    PowerSE{1} = power_dif.' * power_dif;
end % end DAS



%% ------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Detected_powers, Distance, ~, normal]=fun_IAAEst(y_noisy, A, modulus_hat_das,DOAscan,DOA);
if ~normal
%     warning('IAA Abnormal!');
    SqrErr{2} = NaN;
    PowerSE{2} = NaN;
else
    SqrErr{2} = Distance * Distance';
    power_dif = Detected_powers - PowerDOA;
    PowerSE{2} = power_dif.' * power_dif;
end % end IAA








%% --- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Detected_powers, Distance, ~, normal, ~] = fun_SAM1Res(y_noisy, A, modulus_hat_das,DOAscan,DOA);
if ~normal
%     warning('SAMV1 Abnormal!');
    SqrErr{3} = NaN;
    PowerSE{3} = NaN;
else
    SqrErr{3} = Distance * Distance';
    power_dif = Detected_powers - PowerDOA;
    PowerSE{3} = power_dif.' * power_dif;
end % end SAMV-1



%% ---- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Detected_powers_sam3, Distance_sam3, ~, normal_sam3, NoisePower_sam3] = fun_SAM3Res(y_noisy, A, modulus_hat_das,DOAscan,DOA);
if ~normal_sam3
%     warning('SAMV3 Abnormal!');
    SqrErr{4} = NaN;
    PowerSE{4} = NaN;
else
    SqrErr{4} = Distance_sam3 * Distance_sam3';
    power_dif = Detected_powers_sam3 - PowerDOA;
    PowerSE{4} = power_dif.' * power_dif;
    
end % end SAMV-3


%% ---- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Detected_powers, Distance, ~, normal, ~] = fun_SPICEplusRes(y_noisy, A, modulus_hat_das,DOAscan,DOA);
if ~normal
%     warning('SPICEplus Abnormal!');
    SqrErr{5} = NaN;
    PowerSE{5} = NaN;
else
    SqrErr{5} = Distance * Distance';
    power_dif = Detected_powers - PowerDOA;
    PowerSE{5} = power_dif.' * power_dif;
end % end spice+





%% ---- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% provide initial values for the R3_ML and SAMk_ML series...
normal_in = normal_sam3;
initDOA = DOA + Distance_sam3;
initNoisePower = NoisePower_sam3;
initPower = Detected_powers_sam3;




%% Res1_ML
[p_k_mat_r1ML,  ~, ~, Distance_r1ML, normal_r1ML]=...
    fun_Res1_MLRes(y_noisy, A, DOA, initPower, initDOA, initNoisePower, normal_in);
if ~normal_r1ML
    SqrErr{6} = NaN;
    PowerSE{6} = NaN;
else
    SqrErr{6} = Distance_r1ML * Distance_r1ML';
    power_dif = p_k_mat_r1ML.' - PowerDOA;
    PowerSE{6} = power_dif.' * power_dif; % end Res1_ML  
end

%% SAMV-k-ML:
[p_k_mat,  ~, ~, Distance, normal]=...
        fun_SAMk_MLRes(y_noisy, A, DOA, initPower, initDOA, initNoisePower, normal_in, 1); 
if ~normal
    SqrErr{7} = NaN;
    PowerSE{7} = NaN;
else
    SqrErr{7} = Distance * Distance';
    power_dif = p_k_mat.' - PowerDOA;
    PowerSE{7} = power_dif.' * power_dif; %  
end


[p_k_mat,  ~, ~, Distance, normal]=...
        fun_SAMk_MLRes(y_noisy, A, DOA, initPower, initDOA, initNoisePower, normal_in, 3);
if ~normal
    SqrErr{8} = NaN;
    PowerSE{8} = NaN;
else
    SqrErr{8} = Distance * Distance';
    power_dif = p_k_mat.' - PowerDOA;
    PowerSE{8} = power_dif.' * power_dif; %  
end






%% ---- Note: MUSIC must comes last here !!!  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, theta_k_select_music, ~, normal, ~] = ...    
    fun_MUSIC_GridFree_Res(y_noisy, A, modulus_hat_das,DOAscan,DOA);

if ~normal
    SqrErr{9} = NaN;
    PowerSE{9} = NaN;
else
    Distance = theta_k_select_music - DOA;
    SqrErr{9} = Distance * Distance';
    PowerSE{9} = NaN;
end % end MUSIC











% disp('==== end debug function Compute_algos_StdErr.m....');

end % end the function 













