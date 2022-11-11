%% =============== SAM-0 estimates ==========================
function [p_vec, sigma]=SAMV0(Y,A,DAS_init)
% output list: 
% Distance 1 x # source, row vector
% p_vec: # scan point x 1, col vector
% normal: tag, 
% if normal == 1, detecion is Okay
% otherwise normal ==0, detectio failed
%
% Input list: 
% Y: measured data, each col. is one snapshot
% A: steering vector matrix
% DAS_init: initial coefficients estimates by DAS
% Aug 21, 2011 QL
% ---------------------------------------------------

threshold=1e-6;
maxIter= 30;
% colorSet={'r-', 'b-', 'r-.', 'b-.', 'r--', 'b-.', 'r:', 'b:'};

[M, thetaNum]=size(A);
t_samples = size(Y, 2);
R_N = (Y*Y')/t_samples;
sigma =  mean(abs(Y(:)).^2); %1e-6;
% disp(['Actual init sigma in SAMV-0 ==' num2str(sigma) ]);
p_vec_Old = abs(DAS_init).^2; % 
% figure
% plot(20*log10(real(p_vec_Old)))
% hold on
for iterIdx = 1:maxIter
    R = real(A*spdiags(p_vec_Old, 0, thetaNum, thetaNum)*A' + sigma*eye(M));%

    % ----- a little faster 
    Rinv=inv(R);
%     Rinv_RN_Rinv_A = Rinv * R_N * Rinv * A;
%     faktor = real(sum(conj(A).*Rinv_RN_Rinv_A, 1).');
%     p_vec = p_vec_Old.^2 .* faktor;
    sigma = real(trace(Rinv*Rinv*R_N))/real(trace(Rinv*Rinv));
    
    faktor_alternativ = real(diag((A'/R)*R_N*(R\A)));
    p_vec = p_vec_Old.^2 .* faktor_alternativ; % AH recommend this formula instead...

%     plot(20*log10(p_vec))

    if norm(p_vec_Old) == 0 || norm(p_vec_Old-p_vec)/norm(p_vec_Old) < threshold 
%         disp(['=== SAM-0 convgs at iteration ' num2str(iterIdx)]);
        break;
    end
    
    p_vec_Old = p_vec;
end

p_vec = real(p_vec);