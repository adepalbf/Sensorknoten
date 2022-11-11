%% =============== SAM-0 estimates ==========================
function [Detected_powers, Distance, normal] = SAMV0_result(p_vec, DOAscan, DOA)
% output list: 
% Distance 1 x # source, row vector
% normal: tag, 
% if normal == 1, detecion is okay
% otherwise normal ==0, detection failed
%
% Input list: 
% p_vec: # scan point x 1, col vector
% DOAscan: grid
% DOA: truth
% Aug 21, 2011 QL
% ---------------------------------------------------

Numsources =length(DOA);

[~, index]=findpeaks(p_vec, 'sortstr', 'descend');

if length(index) < Numsources
%     warning('Not all peaks detected');
    normal = 0;
    Distance = NaN;
    Detected_powers = NaN;
    return;
end

% ------------ Check whether the detection is right -----
Detected_DOAs = DOAscan(index(1:Numsources));

[Detected_DOAs, IXsort] = sort(Detected_DOAs, 'ascend');
Distance = Detected_DOAs - DOA;
if max(abs(Distance)) > 10
    normal = 0;
    Distance = NaN;
    Detected_powers = NaN;
else
    normal = 1; % detection okay
    % the powers from large value to small value
    Detected_powers = p_vec(index(1:Numsources));
    % sort the power according to the DOA 
    Detected_powers =  Detected_powers(IXsort);
end

end
