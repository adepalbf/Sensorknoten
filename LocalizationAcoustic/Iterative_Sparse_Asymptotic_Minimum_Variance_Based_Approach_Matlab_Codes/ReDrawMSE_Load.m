% script file to redraw the MSE-CRB curves for the SAMV paper draft
% Sep. 8, 2011 by QL
% Based on the Parfor_MC_SAMS.m file
% Load saved MAT Files to replot the figures
%



% ------------------- start input please ...---
cohr_flag = 0; % 1;
M_value = 12;
snap_value
% --------------- end input ---------------




% SNRSetdB
if ~exist('SNRSetdB', 'var')
    disp('Loading falied!... exiting...');
    return;
end

algo_list = {'PER', 'IAA', 'SAMV1', 'SAMV2', 'SPICE+', 'AMV-SML', 'SAMV1-SML', 'SAMV2-SML', 'MUSIC'};
% algo_list = {'PER', 'IAA'};
Num_algos = length(algo_list); % 


SE_mean = zeros(length(SNRSetdB), Num_algos);
Failing_rate = zeros(length(SNRSetdB), Num_algos);

CRBlist = zeros(1, length(SNRSetdB)); % row vec



