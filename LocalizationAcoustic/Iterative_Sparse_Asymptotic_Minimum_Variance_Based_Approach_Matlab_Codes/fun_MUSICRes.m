%% =============== MUSIC estimates ==========================
% note the p in MUSIC is not acutal powers
function [UselessPowerValue, Distance, p, normal, UselessNoisePower]=fun_MUSICRes(y_noisy,A,DAS_init,DOAscan,DOA)
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
% updated: Aug. 21, 2011

UselessPowerValue = NaN; % MUSIC is unable to produce power estimates.
UselessNoisePower = NaN; % not interested in MUSIC noise power estimates.ds
% ---------------------------------------------------
source_no = length(DOA);
const.t_samples = size(y_noisy, 2);
R_hat = y_noisy * y_noisy' /const.t_samples;

[V, D] = eig(R_hat);
d = diag(D);
[dec_d, dec_IX] = sort(d, 'descend');
dec_V = V(:, dec_IX);

% source_no is assumed to be known
G_hat = dec_V(:, source_no+1:end);
MUSIC_out = 1./ (sum(conj(A) .* (G_hat * G_hat' * A), 1));
% convert to col vec, and discard small imag part
p = real(MUSIC_out.'); 


[pks index]=findpeaks(p, 'sortstr', 'descend');

if length(index) < source_no
%     warning('Not all peaks detected');
    normal = 0;
    Distance = NaN;
%     p = NaN;
    return;
end

% ------------ Check whether the detection is right -----
Detected_DOAs = DOAscan(index(1:source_no));
Detected_DOAs = sort(Detected_DOAs, 'ascend');
Distance = Detected_DOAs - DOA;
if max(abs(Distance)) > 10
%     warning('Failed Detection by MUSIC, this simulation is abnormal!');
    normal = 0;
    Distance = NaN;
%     p = NaN;
else
    normal = 1; % detection okay
end

end

