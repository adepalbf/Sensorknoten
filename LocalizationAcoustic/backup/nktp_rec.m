function x=nktp_rec(SAMPLES, FREQ)
% x=nktp_record(SAMPLES, FREQ)
% get real data from the four microphone arrays.
% inputs:
%  1- SAMPLES: number of samples to acquire.
%  2- FREQ: the sampling frequency.
% output:
%  x: the output signal at the microphones.

recObj = audiorecorder(FREQ,8,2,1);
disp('Recording...')
recordblocking(recObj, SAMPLES/FREQ);
disp('...finished')
x(:, 1:2) = getaudiodata(recObj);

% AI1 = analoginput( 'winsound',1);
% AI2 = analoginput( 'winsound',2);
% AI3 = analoginput( 'winsound',3);
% AI4 = analoginput( 'winsound',4);
% ch1 = addchannel(AI1,1:2);
% ch2 = addchannel(AI2,1:2);
% ch3 = addchannel(AI3,1:2);
% ch4 = addchannel(AI4,1:2);
% set([AI1 AI2 AI3 AI4], 'SampleRate',FREQ)
% set([AI1 AI2 AI3 AI4], 'SamplesPerTrigger',SAMPLES)
% start([AI1 AI2 AI3 AI4]) %start and trigger analog inputs
% wait([AI1 AI2 AI3 AI4], SAMPLES/FREQ) %wait for analog inputs
% x(:,1:2)=getdata(AI1); %get recorded data
% x(:,3:4)=getdata(AI2);
% x(:,5:6)=getdata(AI3);
% x(:,7:8)=getdata(AI4);
% 
% delete([AI1 AI2 AI3 AI4]);
% clear AI1 AI2 AI3 AI4;

end