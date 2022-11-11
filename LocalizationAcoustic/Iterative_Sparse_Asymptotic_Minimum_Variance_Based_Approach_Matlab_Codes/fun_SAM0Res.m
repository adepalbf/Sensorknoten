%% =============== SAM-0 estimates ==========================
function [Detected_powers, Distance, p_vec, normal, noisepower]=fun_SAM0Res(Y,A,DAS_init,DOAscan,DOA)
% Eliminate the for-loop 
%
% simply put the private function out
% Updated Jun 19, 2011 QL
% ---------------------------------------------------
% output list: 
% Distance 1 x # source, row vector
% p_vec: # scan point x 1, col vector
% normal: tag, 
% if normal == 1, detecion is Okay
% otherwise normal ==0, detectio failed
%
%
% Input list: 
% Y: measured data, each col. is one snapshot
% A: steering vector matrix
% DAS_init: intitial coefficients estimates by DAS
% DOAscan: grid
% DOA: truth
% Aug 21, 2011 QL
% ---------------------------------------------------

Numsources =length(DOA);
threshold=1e-6;
maxIter= 30;
% colorSet={'r-', 'b-', 'r-.', 'b-.', 'r--', 'b-.', 'r:', 'b:'};

[M thetaNum]=size(A);
t_samples = size(Y, 2);
R_N = (Y*Y')/t_samples;
sigma =  mean(abs(Y(:)).^2); %1e-6;
% disp(['Actual init sigma in SAMV-0 ==' num2str(sigma) ]);
p_vec_Old = abs(DAS_init).^2; % 


for iterIdx = 1:maxIter
    R =  A*spdiags(p_vec_Old, 0, thetaNum, thetaNum)*A' + sigma*eye(M);% 
    Rinv=inv(R);
%     p_vec = zeros(thetaNum, 1);
%     for snapIdx = 1:t_samples
%         y = Y(:, snapIdx);
%         Rinv_y = Rinv*y;
% %         Rinv_A = Rinv*A;
% %         diag_A_Rinv_A = sum(conj(A).*Rinv_A, 1).';
%         coeff = A'*Rinv_y;%./real(diag_A_Rinv_A);
%         p_vec = p_vec + p_vec_Old.^2.*abs(coeff).^2;
%     end
    
    
    % ----- a little faster 
    Rinv_RN_Rinv_A = Rinv * R_N * Rinv * A;
    p_vec = p_vec_Old.^2 .* sum(conj(A).*Rinv_RN_Rinv_A, 1).';
    
%     p_vec = p_vec_Old.^2 .* (diag(A'*Rinv*R_N*Rinv*A)); % AH recommend this formula instead...

%     % ----- this is slow
%     Rinv_RN_Rinv = Rinv * R_N * Rinv;
%     p_vec = p_vec_Old.^2 .* diag(A' * Rinv_RN_Rinv * A );
    
    
    
%     p_vec = p_vec/t_samples;
    
    sigma = real(trace(Rinv*Rinv*R_N))/real(trace(Rinv*Rinv));
    
    if norm(p_vec_Old-p_vec)/norm(p_vec_Old)<threshold
%         disp(['=== SAM-0 convgs at iteration ' num2str(iterIdx)]);
        break;
    end
    
    p_vec_Old = p_vec;
end

p_vec = real(p_vec);
[pks index]=findpeaks(p_vec, 'sortstr', 'descend');

if length(index) < Numsources
%     warning('Not all peaks detected');
    normal = 0;
    Distance = NaN;
    % keep the p_vec values for outbound plotting...
%     p_vec = NaN;
    Detected_powers = NaN;
    noisepower = sigma;
    return;
end


% ------------ Check whether the detection is right -----
Detected_DOAs = DOAscan(index(1:Numsources));

[Detected_DOAs, IXsort] = sort(Detected_DOAs, 'ascend');
Distance = Detected_DOAs - DOA;
if max(abs(Distance)) > 10
%     warning('Failed Detection by SAM-0, this simulation is abnormal!');
    normal = 0;
    Distance = NaN;
    % keep the p_vec values for outbond plotting figures
%     p_vec = NaN;
    Detected_powers = NaN;
else
    normal = 1; % detection okay
    % the powers from large value to small value
    Detected_powers = p_vec(index(1:Numsources));
    % sort the power according to the DOA 
    Detected_powers =  Detected_powers(IXsort);
end

noisepower = sigma;


end
