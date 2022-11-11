%% =============== SPICE estimates based on SPICE_LIKES NEW codes ========================== %
function [Detected_powers, Distance, p, normal, noisepower]=fun_SPICENew(Y,A,DAS_init,DOAscan,DOA)
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
% Sep. 5, 2011 by QL based on official SPCIE+ codes...
% ---------------------------------------------------

Numsources =length(DOA);

maxIter=25;
thetaNum = size(A, 2);




[~,p] = SPICE_MS(A,Y,maxIter); % note here is not SPICE+, p is of length K+M

p = real(p(1:thetaNum));

[pks index]=findpeaks(p, 'sortstr', 'descend');

if length(index) < Numsources
%     warning('Not all peaks detected');
    normal = 0;
    Distance = NaN;
%     p = NaN;
    Detected_powers = NaN;
    noisepower = sigma;
    return;
end

% ------------ Check whether the detection is right -----
Detected_DOAs = DOAscan(index(1:Numsources));

[Detected_DOAs, IXsort] = sort(Detected_DOAs, 'ascend');
Distance = Detected_DOAs - DOA;
if max(abs(Distance)) > 10
%     warning('Failed Detection by SPICE, this simulation is abnormal!');
    normal = 0;
    Distance = NaN;
%     p = NaN;
    Detected_powers = NaN;
else
    normal = 1; % detection okay
    % the powers from large value to small value
    Detected_powers = p(index(1:Numsources));
    % sort the power according to the DOA 
    Detected_powers =  Detected_powers(IXsort);
end

noisepower = NaN; %sigma;



end
