%% =============== IAA estimates ==========================
function [Detected_powers, Distance, p_vec, normal]=fun_IAAEst(Y,A,DAS_init,DOAscan,DOA)
% Put the private function out 
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
% % Verified June 8, 2011 QL
% ---------------------------------------------------

Numsources =length(DOA);
threshold=1e-6;
maxIter=45;
% colorSet={'r-', 'b-', 'r-.', 'b-.', 'r--', 'b-.', 'r:', 'b:'};

[M thetaNum]=size(A);
t_samples = size(Y, 2);
% RE = 1/t_samples*(Y*Y');
% sigma = 1e-10;
p_vec_Old = abs(DAS_init).^2; % 
% coeffOld = DAS_init;

for iterIdx = 1:maxIter
    R =  A*spdiags(p_vec_Old, 0, thetaNum, thetaNum)*A';% + sigma*eye(M);% R= R+mean(diag(R))*1e-12*eye(size(R,1));
    Rinv=inv(R);
    p_vec = zeros(thetaNum, 1);
    for snapIdx = 1:t_samples
        y = Y(:, snapIdx);
        Rinv_y = Rinv*y;
        Rinv_A = Rinv*A;
        diag_A_Rinv_A = sum(conj(A).*Rinv_A, 1).';
        coeff = A'*Rinv_y./real(diag_A_Rinv_A);
        p_vec = p_vec + abs(coeff).^2;
%         sigma = real(trace(Rinv*Rinv*RE))/real(trace(Rinv*Rinv));
    end
    p_vec = p_vec/t_samples;
    
    if norm(p_vec_Old-p_vec)/norm(p_vec_Old)<threshold
%         disp(['IAA convgs at iteration ' num2str(iterIdx)]);
        break;
    end
    p_vec_Old = p_vec;
end
[~, index]=findpeaks(p_vec, 'sortstr', 'descend');

if length(index) < Numsources
%     warning('Not all peaks detected');
    normal = 0;
    Distance = NaN;
%     p_vec = NaN;
    Detected_powers = NaN;
    return;
end

% ------------ Check whether the detection is right -----
Detected_DOAs = DOAscan(index(1:Numsources));

[Detected_DOAs, IXsort] = sort(Detected_DOAs, 'ascend');
Distance = Detected_DOAs - DOA;
% if max(abs(Distance)) > 10
% %     warning('Failed Detection by IAA, this simulation is abnormal!');
%     normal = 0;
%     Distance = NaN;
% %     p_vec = NaN;
%     Detected_powers = NaN;
% else
%     normal = 1; % detection okay
%     % the powers from large value to small value
%     Detected_powers = p_vec(index(1:Numsources));
%     % sort the power according to the DOA 
%     Detected_powers =  Detected_powers(IXsort);
% end


    normal = 1; % detection okay
    % the powers from large value to small value
    Detected_powers = p_vec(index(1:Numsources));
    % sort the power according to the DOA 
    Detected_powers =  Detected_powers(IXsort);


end
