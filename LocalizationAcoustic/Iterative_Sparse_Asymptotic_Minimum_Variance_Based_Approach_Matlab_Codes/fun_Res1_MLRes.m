%% =============== Refined: Result 1 {eq. 4, 5} + fminsearch init by SAMV-3 ==========================
function [p_k_mat,  noisePower, theta_k_mat, Distance, normal]=...
    fun_Res1_MLRes(Y,A,DOA, initPower, initDOA, initNoisePower, normal_in) % OK...
% Assume: # of sources are known
% a lot of debugging options, fminsearch smaller TolX, TolFun: 1e-6 
% ---------------------------------------------------
% output list: 
% p_k_mat: the estimated power of the K sources
% theta_k_mat: the estimated source locations, grid-free methods
% normal: always 1
% noisePower: refined noise power update by the remark 3 algos...
%
% Input list: 
% Y: measured data, each col. is one snapshot
% A: steering vector matrix
% DOA: truth DOA values
% initPower: initial signal power, provided by SAMV-3
% initDOA: initial DOA estiamtes, provided by SAMV-3
% initNoisePower: initial noise estimates, provided by SAMV-3
% normal_in: flag whether the intial SAMV-3 is normal
%
%
%  Aug. 21, 2011 by QL
% ---------------------------------------------------
if ~normal_in
    normal = 0; % initialization failed!
    p_k_mat = NaN;
    noisePower = NaN;
    theta_k_mat = NaN;
    Distance = NaN;
    return;
end
    


%% === 
threshold = 1e-6;
Numsources =length(DOA);
t_samples = size(Y, 2);
[M thetaNum]=size(A);

R_N = (Y*Y')/t_samples;

% %% ===== begin inside initialization.... 
% modulus_hat_das  = sum(abs(A'*Y /M), 2 )/t_samples;
% [Detected_powers_current, Distance_current, ~, normal_init, NoisePower_init] = fun_SAM3Res(Y,A,modulus_hat_das,DOAscan,DOA);
% % it may fail, so check the normal values... here!
% if ~normal_init % intial SAMV-3 failed....
%     p_k_mat=NaN;
%     sigma_updata=NaN; 
%     theta_k_mat=NaN; 
%     Distance=NaN; 
%     normal=0; 
%     noisePower=NaN;
%     return;
% end
% detected_DOAs_current = DOA + Distance_current; % intial DOA given by SAM3
% sigma_current = NoisePower_init; 
% %% ==== end inside initialization.....

% use parameter to do initialization

Detected_powers_current = initPower;
detected_DOAs_current = initDOA;
sigma_current = initNoisePower;

%% --- we dynamically update the steering vector now, not using the DOAscanning stv matrix A
Dist = ones(1, M-1); 
DistTmp = cumsum([0 Dist]);

% ----
maxIter=60;
for iterIdx = 1:maxIter
    A_current = exp(1j*pi*DistTmp' * cos(detected_DOAs_current *pi/180) );
    
    R = A_current * diag(Detected_powers_current) * A_current' + sigma_current * eye(M);
    Rinv = inv(R);
    
    sigma_updata = real(trace(Rinv*Rinv*R_N) + sigma_current * trace(Rinv*Rinv) - trace(Rinv) )/real(trace(Rinv*Rinv));
%         if sigma_updata < 0
%             disp('Sigma < 0 and setting it to zero NOW... ');
%             sigma_updata = 0;
%         end
    
    % store p_k, theta_k for updating
    p_k_mat = zeros(1, Numsources);
    theta_k_mat = zeros(1, Numsources);
    for source_ind = 1: Numsources
        Q_k = R - Detected_powers_current(source_ind) * ...
            A_current(:, source_ind) * A_current(:, source_ind)';
        p_k = A_current(:, source_ind)' * Rinv * R_N * Rinv * A_current(:, source_ind)/...
            ((A_current(:, source_ind)' * Rinv * A_current(:, source_ind)  )^2) + ...
            Detected_powers_current(source_ind) - 1/( A_current(:, source_ind)' * Rinv* A_current(:, source_ind)  );
        
        p_k = real(p_k); % discard the tiny imaginary part
        
%         if p_k < 0 
%             disp('p_k < 0 and setting it to zero NOW...');
%             p_k = 0;
%         end
        
        % save p_k values
        p_k_mat(source_ind) = p_k;
        

        % Search 1D using fminsearch function...

        [theta_k_update ] = fminsearch(@(theta) myfun(theta,p_k, Q_k, R_N, M),   ...
            detected_DOAs_current(source_ind), optimset('TolX', 1e-6, 'TolFun', 1e-6)    ); 

         theta_k_mat(source_ind) = theta_k_update;
    end
    
    

    
    % exit the loop using termination criterions
    % threshold
    if norm(detected_DOAs_current - theta_k_mat) < threshold % practical convg
%         disp(['Thetas Convg @ iter ' num2str(iterIdx)  ]);
        break;
    end
    
    
    % updating formula
    Detected_powers_current = p_k_mat;
    sigma_current = sigma_updata;
    detected_DOAs_current = theta_k_mat;
    

end % end iters...



% after the iterations, the useful information is stored in:
% p_k_mat, sigma_updata, theta_k_mat these variables...


Distance = theta_k_mat - DOA;
noisePower = sigma_updata;

% % check it is not too far away...
% if max(abs(Distance)) > 10 % deg
%     normal = 0;
%     disp('Res1_ML failed due to Distance > 10! ');
% else
%     normal = 1;
% end

normal = 1;



end







%% searching function used by fminsearch....
function f = myfun(theta, p_k, Q_k, R_N, M)
Dist = ones(1, M-1); 
DistTmp = cumsum([0 Dist]);
a_k = exp(1j*pi*DistTmp' * cos(theta *pi/180) );
b_k = Q_k\a_k;
beta_k = 1/ (1 + p_k * a_k' * b_k  );
f = log(1 + p_k * a_k' * b_k )  -  p_k * beta_k * (b_k' * R_N * b_k); 
f = real(f);
end





