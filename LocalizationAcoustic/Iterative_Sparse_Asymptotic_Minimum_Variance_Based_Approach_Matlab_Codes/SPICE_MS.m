% this is the offical codes of mulitple-snapshot SPICE in the to-appear
% (Not published) SPICE-and-LIKES paper
%
% This codes is obtained at the begining of Aug. 2011, possibly written by
% Prabhu Babu, Uppsala Univ. Sweden.
function [beta,p] = SPICE_MS(A,Y,iterations)

s    = size(Y,2);
N    = size(Y,1);
Rh   = Y*Y'/s;
B    = [A,eye(N)];
p    = diag(abs(B'*Rh*B));
R    = B*diag(p)*B';
Rinv = inv(R);
w    = sqrt(sum(abs(B).^2,1));
w    = w';

for iter = 1:iterations
    beta       = diag(p)*B'*inv(R)*Y;
    pprev      = p;
    p          = sqrt(sum(abs(beta).^2,2))./w;
    R          = B*diag(p)*B';
    if(norm(pprev-p)/norm(p) < 10^-3)
        break;
    end
end