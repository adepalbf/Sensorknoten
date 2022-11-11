
%The antenna operating frequency is 150 MHz.

fc = 150.0e6;            
phased_array = phased.ULA('NumElements',4,'ElementSpacing',0.15);
%Create the arriving signals at the ULA. The true direction of arrival of the first signal is 10° in azimuth and 20° in elevation.
%The direction of the second signal is 60° in azimuth and -5° in elevation.
 
fs = 8000.0;                                 
t = (0:1/fs:1).';       
sig1 = cos(2*pi*t*300.0);
sig2 = cos(2*pi*t*400.0);
sig = collectPlaneWave(phased_array,[sig1 sig2],[10 20; 60 -5]',fc);
noise = 0.1*(randn(size(sig)) + 1i*randn(size(sig)));

%Estimate the DOAs. 

estimator = phased.MUSICEstimator('SensorArray',phased_array,... 
    'OperatingFrequency',fc,...
    'DOAOutputPort',true,'NumSignalsSource','Property',...
    'NumSignals',2);
[y,doas] = estimator(sig + noise);
doas = broadside2az(sort(doas),[20 -5])

%Plot the MUSIC spectrum.

plotSpectrum(estimator,'NormalizeResponse',true)