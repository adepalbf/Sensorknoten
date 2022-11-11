function [beta, Detected_powers]  = AOA_MUSIC(Y, ULA)
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
        [beta(ii), spec, specang] = musicdoa(covmat,N, 'ElementSpacing',d);
    
    plot(specang,10*log10(spec))
    xlabel('Arrival Angle (deg)')
    ylabel('Magnitude (dB)')
    title('MUSIC Spectrum')
    grid
    
    hold all
    [Detected_powers(ii), index] = max(spec);%rxy
    catch
    end
end

