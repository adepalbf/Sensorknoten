%% RELAX extension for the simple 1D SAMs algos
function [NewPs, NewDOAs] = PlusRelax(DOAinit, beta_init, y_noisy   ) 

source_no = length(DOAinit);  % must be two here
MaxRELAX_no = 45;
%% ongoing here @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

y_temp    = zeros(size(y_noisy));
for repeat = 1:MaxRELAX_no
    for kk = 1: source_no
        for Find = 1: Fbins
            sigset = sigsetFB(Find, :); % row vector here
            sigset_temp = sigset.'; % col vector
            
            sigset_temp(kk) = 0; % only useful for 2 sources
            
            % generate stv
            deltaR = funCalDeltaR(sonar, rangeSet, (DOAset*pi/180), elevAngleSet, algParam.modelType);
            stv = exp(1j*2*pi*freqSet(Find)*deltaR/const.c);
            
            y_temp(:, Find) = y_noisy(:, Find) - stv * sigset_temp;
        end % end f-bins....

       % on-going 
       [theta, val] = ...
           fminsearch(@(theta) myfun(theta, y_temp, Fbins, sonar, algParam), DOAset(kk), optimset('TolX', 2e-4  ));

       DOAsetnew(kk) = theta;
       
%        % ==== for debug only ---- @@@@@@@@@@@@@
%        disp(['theta new values after search : ' num2str(theta) ]);
%        % ----------------------@@@@@@@@@@@@@@
       
       for Find = 1: Fbins
           % update coeff's...
           % ---- generate updated stv 
           deltaR = funCalDeltaR(sonar, rangeSet, (theta *pi/180), elevAngleSet, algParam.modelType);
           stv_temp = exp(1j*2*pi*freqSet(Find)*deltaR/const.c);
           sigsetnewFB(Find, kk) = stv_temp' * y_temp(:, Find) / sonar.M;
       end % end Fbins...
           

    end % end kk --- source_no
       sigsetFB = sigsetnewFB;
       
       % early termination 
       if sum(abs(DOAset - DOAsetnew))< 1e-5 % not changed!!
%        if norm(DOAset - DOAsetnew) / norm(DOAset) < 1e-6 % not changed!!
           disp(['WB-RELAX convgs after repeat = ' num2str(repeat) ]);
           break;
       end
       DOAset = DOAsetnew;
end % end repeat

% output parameters....
NewDOAs = DOAset;

% sigsetFB is 12 x 2 matrix...
NewPs =   diag(sigsetFB' * sigsetFB).' / sonar.sampNumPerBlk; % row vecotor

