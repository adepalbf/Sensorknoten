filename = 'C:\Users\Adepa\Documents\MCUXpressoIDE_11.5.0_7232\workspace\sensor_node_git\fc1_buffer_7.csv';
nByte = 4;
nbit = 8*nByte;

nBuffer = 960 /4*nByte;
q = quantizer('fixed', 'nearest', 'wrap', [nbit 0]);% quantizer object for num2hex function  
FID = fopen(filename);
dataFromfile = textscan(FID, '%s');% %s for reading string values (hexadecimal numbers)
dataFromfile = dataFromfile{1};
data2 = cell(length(dataFromfile)/4*nByte, 1);
for ii = 1:length(dataFromfile)/4
    data2((ii-1)*nByte+(1:nByte)) = dataFromfile((ii-1)*4+(nByte:-1:1));
end
dataFromfile = data2;
left_channel = cell(fix(nBuffer/2/nByte)-1,1);
right_channel = cell(fix(nBuffer/2/nByte)-1,1);
offset = 0;
for ii = 1:fix(nBuffer/2/nByte)-1
    left_channel{ii} = [dataFromfile{offset+(1:nByte)+(ii-1)*2*nByte}];
    right_channel{ii} = [dataFromfile{offset+nByte+(1:nByte)+(ii-1)*2*nByte}];
end
leftData = hex2num(q, left_channel);
rightData = hex2num(q, right_channel);

% dataDec   = base2dec((left_channel), 24);
% index     = dataDec >= 2^23;
% dataDec(index) = dataDec(index) - 2^24;


leftData = cell2mat(leftData);
rightData = cell2mat(rightData);

% leftData = dataDec;
%%

data = [leftData,rightData];
jump_idx = 90;
data = [data(jump_idx:end,:); data(1:jump_idx-1,:)];
figure
plot(data)
% plot(diff([leftData,rightData]))
%%
figure
pwelch(data, [],[],[],48000);
% decData = hex2num(q, dataFromfile);
% decData = cell2mat(decData);

%%
fclose(FID);