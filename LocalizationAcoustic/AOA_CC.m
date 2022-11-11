
function [beta, Detected_powers]  = AOA_CC(Y, dn_max)
% calculate the generalized cross-correlation

M = size(Y,2);
if isscalar(dn_max)
    dn_max = repmat(dn_max, M-1, 1);
end
Detected_powers = zeros(M-1, 1);
beta = zeros(M-1, 1);
figure
 
for ii = 1:M-1
    max_lag = dn_max(ii);
    [rxy, lags] = xcorr(Y(:,1), Y(:,ii+1),max_lag*2, 'unbiased');
    lags_interp = lags(1):0.01:lags(end);
    rxy_interp = interp1(lags, rxy, lags_interp, 'spline');
%     stem(lags, rxy)
%     hold on
    plot(lags_interp, rxy_interp, 'DisplayName', "r_{"+1+""+(ii+1)+"}")
    hold on
    
    [Detected_powers(ii), index] = max(rxy_interp);
    
    beta(ii) = acos(lags_interp(index)/max_lag);
    stem(lags_interp(index), rxy_interp(index), 'DisplayName', "Local AOA: "+round(180/pi*beta(ii),2))
%     stem(max_lag, rxy(index))
    
    % calculate the angles
end
grid on
xlabel("Lags")
ylabel("Correlation coefficient")
legend('Location', 'best')
