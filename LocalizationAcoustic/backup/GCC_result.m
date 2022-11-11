%% =============== SAM-0 estimates ==========================
function [Distance, normal] = GCC_result(Detected_DOAs, DOA)
% output list:
% Distance 1 x # source, row vector
% normal: tag,
% if normal == 1, detecion is okay
% otherwise normal ==0, detection failed
%
% Input list:
% Detected_DOAs
% DOA: truth
% Aug 21, 2011 QL
% ---------------------------------------------------

Numsources =length(DOA);
N_sources = min(Numsources, length(Detected_DOAs));

Distance = NaN(N_sources);

for ii = 1:N_sources
    [Distance(ii), idx] = min(abs(Detected_DOAs - DOA(ii)));
    Detected_DOAs(idx) = [];
end
if max(abs(Distance)) > 10
    disp('High error')
    normal = 0;
else
    normal = 1; % detection okay
end

end
