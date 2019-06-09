clc; clear all;

channel = {'Multipath','Fading','ByPass', 'AWGN'};

% Parameters
P.NumberOfFrames   = 50;
P.NumberOfBits     = 172;
P.Q_Ind = 12;

% Only applies for fading and multipath
lengthch=[3,4,5,6,10];
P.ChannelLength = 3; 

P.RakeFingers = 3; 
if (P.RakeFingers > P.ChannelLength)
    error('Fingers has to be smaller or equal to channels !')
end

% Only applies for fading
P.CoherenceTime = 19200/1000; % A thousandth of a second 

P.HadLen = 64; % Length of Hadamard Sequence

P.K = 9; % Length of convolutional encoder
P.ConvSeq = [753 561]; % Rate 1/2
P.Rate = length(P.ConvSeq);

P.LongCodeLength = 42; % PN Sequence

P.SequenceMask = [1,1,0,0,0,1,1,0,0,0, randi([0 1],1,32)];

P.SNRRange = -50:5:0; % SNR Range to simulate in dB

%changing channels
figure;
hold on;
grid minor;
for i=1:4
P.ChannelType=channel{i} 

BER = SISOsimulator(P);

if strcmp(P.ChannelType, 'Multipath') | strcmp(P.ChannelType, 'Fading')
    simlab = sprintf('%s - Paths: %d - Fingers : %d' ,P.ChannelType,P.ChannelLength, P.RakeFingers);
else
    simlab = sprintf('%s' ,P.ChannelType);
end

%figure;
%hold on;
semilogy(P.SNRRange,BER,'-','DisplayName',simlab)

xlabel('SNR','FontSize',12,'FontWeight','bold');
ylabel('BER','FontSize',12,'FontWeight','bold');
xlim([min(P.SNRRange) max(P.SNRRange)]);

legend('-DynamicLegend');
end

%changing channel length
figure;
hold on;
P.ChannelType='Multipath'

for j=1:length(lengthch)
    P.ChannelLength=lengthch(j);
    
    BER = SISOsimulator(P);

if strcmp(P.ChannelType, 'Multipath') | strcmp(P.ChannelType, 'Fading')
    simlab = sprintf('%s - Paths: %d - Fingers : %d' ,P.ChannelType,P.ChannelLength, P.RakeFingers);
else
    simlab = sprintf('%s' ,P.ChannelType);
end

txt = ['Length = ',num2str(lengthch(j))];
semilogy(P.SNRRange,BER,'DisplayName',txt)

xlabel('SNR','FontSize',12,'FontWeight','bold');
ylabel('BER','FontSize',12,'FontWeight','bold');
xlim([min(P.SNRRange) max(P.SNRRange)]);

legend('-DynamicLegend');
    
end

