%% =============== MUSIC -Grid Free Using fminsearch  ==========================
% Also suitable for computing MSE, etc...
%
% Aug. 21, 2011  QL...

% note the p in MUSIC is not acutal powers
function [pseu_p_select, theta_k_select, pseu_p, normal, noisepower] = ...
    fun_MUSIC_GridFree_Res(y_noisy,A,DAS_init,DOAscan,DOA)
% ---------------------------------------------------
% output list: 
%
%
%
%
%
% Input list: 
% Y: measured data, each col. is one snapshot
% A: steering vector matrix
% DAS_init: intitial coefficients estimates by DAS
% DOAscan: grid
% DOA: truth
% updated: Aug. 21, 2011


noisepower = NaN; % not relevant values.
M = size(A, 1);
% ---------------------------------------------------
source_no = length(DOA);
const.t_samples = size(y_noisy, 2);
R_hat = y_noisy * y_noisy' /const.t_samples;

[V, D] = eig(R_hat);
d = diag(D);
[~, dec_IX] = sort(d, 'descend');
dec_V = V(:, dec_IX);

% source_no is assumed to be known
G_hat = dec_V(:, source_no+1:end);
MUSIC_out = 1./ (sum(conj(A) .* (G_hat * G_hat' * A), 1));
% convert to col vec, and discard small imag part
p = real(MUSIC_out.'); 


[~, index]=findpeaks(p, 'sortstr', 'descend');

if length(index) < source_no
%     warning('Not all peaks detected');
    normal = 0;
    pseu_p = p;
    theta_k_select =  NaN;
    pseu_p_select = NaN;
    return;
end

% if proceed, right Num of peaks are detected...

% ------------ Check whether the detection is right -----
Detected_DOAs = DOAscan(index(1:source_no));
Detected_DOAs = sort(Detected_DOAs, 'ascend');
Distance = Detected_DOAs - DOA;
% if max(abs(Distance)) > 30
% %     warning('Failed Detection by MUSIC, this simulation is abnormal!');
%     normal = 0;
%     pseu_p = p;
%     theta_k_select =  NaN;
%     pseu_p_select = NaN;
%     return;
% else
%     normal = 1; % detection okay
%     pseu_p = p;
%     % start grid-free search process...
% end


normal = 1; % detection okay
pseu_p = p;



% begin Grid-free Searching
pseu_p_select = zeros(1, source_no);
theta_k_select = zeros(1, source_no);
for ind = 1: source_no
    [theta_k_select(ind),  pseu_p_select(ind)] = ...
        fminsearch(@(theta)   musicSearch(theta, G_hat, M), Detected_DOAs(ind), ...
        optimset('TolX', 1e-6, 'TolFun', 1e-6)   ); 
end


end


%% ----- 
function f = musicSearch(theta, G_hat, M)
% construct a_k inside the function
Dist = ones(1, M-1); 
DistTmp = cumsum([0 Dist]);
a_k = exp(1j*pi*DistTmp' * cos(theta *pi/180) );
f = 1 / (a_k' * G_hat * G_hat' * a_k);
f = - real(f);
end % end musicSearch function
