function [theta, Detected_powers, rxy, I]  = AOA_GCC_AV(Y, N, dn_max, I_old, alpha)
% calculate the generalized cross-correlation

[rxy, I] = genxcorr_av(Y, N, I_old, alpha);

% find maximum position
% index = find(rxy==max(rxy), 1);
% index = zeros(size(rxy, 2),1);
for ii = 1:size(rxy, 2)
    
    [~, index] = findpeaks(rxy(:,ii), 'sortstr', 'descend');
    
    delay = index;
    % shift delay around center
    delay(delay > N/2) = delay(delay > N/2) - N;
    % remove all outside max. delay
    
%     idx = abs(delay) > dn_max(ii);
%     delay(idx) = [];
%     index(idx) = [];

    % take the strongest remaining
    delay = delay(1);
    index = index(1);
    % the powers from large value to small value
    Detected_powers(ii) = rxy(index,ii);
    
    % calculate the angles
    theta(ii) = asin(delay/dn_max(ii));
    
    % sort the DOA 
%     [theta{ii}, IXsort] = sort(angles, 'ascend');
    % sort the power according to the DOA 
%     Detected_powers{ii} =  Powers(IXsort);
end
