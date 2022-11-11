function [beta, Detected_powers, rxy]  = AOA_GCC(Y, dn_max, N)
% calculate the generalized cross-correlation
gamma = 1e-7;

M = size(Y,2);
if isscalar(dn_max)
    dn_max = repmat(dn_max, M-1, 1);
end
Detected_powers = zeros(M-1, 1);
beta = zeros(M-1, 1);
for ii = 1:M-1
    max_lag = dn_max(ii);
    [rxy, lags] = xcorr(Y(:,ii+1),Y(:,1), N, 'unbiased');
    %% generalized
    Sxy = fft(rxy);
    Rxy = real(ifft(Sxy./(abs(Sxy)+gamma)));
    
    %%
    idx = abs(lags) > max_lag;
    rxy(idx) = [];
    Rxy(idx) = [];
    lags(idx) = [];
    
%     plot(lags, rxy/max(rxy))
%     hold all
%     plot(lags, Rxy/max(Rxy))
    
    [~, index] = max(Rxy);
    Detected_powers(ii) = rxy(index);
    % calculate the angles
    beta(ii) = asin(lags(index)/dn_max(ii));
end

% 
% if isscalar(dn_max)
%     dn_max = repmat(dn_max, size(rxy, 2), 1);
% end
% % find maximum position
% % index = find(rxy==max(rxy), 1);
% % index = zeros(size(rxy, 2),1);
% for ii = 1:size(rxy, 2)
%     
%     [~, index] = findpeaks(rxy(:,ii), 'sortstr', 'descend');
%     
%     delay = index;
%     % shift delay around center
%     delay(delay > N/2) = delay(delay > N/2) - N;
%     % remove all outside max. delay
%     idx = abs(delay) > dn_max(ii);
%     delay(idx) = [];
%     index(idx) = [];
%     % take the strongest remaining
%     delay = delay(1);
%     index = index(1);
%     % the powers from large value to small value
%     Detected_powers(ii) = rxy(index,ii);
%     
%     % calculate the angles
%     theta(ii) = asin(delay/dn_max(ii));
%     
%     % sort the DOA 
% %     [theta{ii}, IXsort] = sort(angles, 'ascend');
%     % sort the power according to the DOA 
% %     Detected_powers{ii} =  Powers(IXsort);
% end
