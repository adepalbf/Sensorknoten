
function [beta, Detected_powers]  = AOA_GCC_PHAT(Y, dn_max)
% calculate the generalized cross-correlation with phase transform

[selected_lag, rxy] = gccphat(Y(:,2:end), Y(:,1));
Detected_powers = max(rxy);

beta = acos(selected_lag(:)./dn_max(:));
