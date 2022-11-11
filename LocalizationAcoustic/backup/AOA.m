function theta = AOA(X1, X2, dn_max, N)
% calculate the cross-correlation
rxy = xycorr(X1, X2, N);
% find maximum position
dn = find(rxy==max(rxy));
% flip results from the high side of the spectrum to the low side
if dn>N/2
    dn = dn-N;
end
% calculate the angle
theta = asin(dn/dn_max);
