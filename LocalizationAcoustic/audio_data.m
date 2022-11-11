function y = audio_data(var_data, fs, varargin)

switch var_data
    case 1
        T = varargin{1};
        f = varargin{2};
        % sine
        t = 1:T*fs;
        y = sin(f*t/fs*2*pi);
    case 2
        T = varargin{1};
        % white Gaussian signal
        y = randn(T*fs,1);
    case 3
        % beepbeep
        [y, Fs] = audioread('beepbeep.wav');  
        y = resample(y,fs,Fs);
end
