
function [beta, Detected_powers]  = AOA_ESPRITE(Y, ULA)
% calculate the generalized cross-correlation

M = size(Y,2);
figure
N = 1;
beta = nan(M-1,1);
Detected_powers = zeros(M-1,1);
for ii = 1:M-1
    d = sqrt(sum((ULA(:,1) - ULA(:,ii+1)).^2));
    covmat = cov(Y(:,1), Y(:,ii+1));
    try
        beta(ii) = espritdoa(covmat,N, 'ElementSpacing', d);
    catch
    end
end