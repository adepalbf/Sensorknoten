function [A_real, A_scan] = scanning_arrays(f, phi_source, mics_xy, phi_scan)

v = 343;


lambda = v/f;
k = 2*pi/lambda*[cos(phi_source); cos(phi_source)]; % sources in columns

% Areal = exp(1i*pi*DistTmp' * cos(DOA*pi/180) ); % real steering vector matrix
A_real = exp(1i*mics_xy*k);

% steering vector matrix w.r.t all possible scanning DOA's
k = 2*pi/lambda*[cos(phi_scan); sin(phi_scan)]; % sources in columns
A_scan = exp(-1i*mics_xy*k);
% A = exp(1i*pi*DistTmp' * cos(DOAscan*pi/180) );
