% Wireless Receivers Project:
% Anael Buchegger, Tim Tuuva, David Sanchez
%
% CDMA Parameter File
% SISO, 1 User only (simplest case)
%
% Telecommunications Circuits Laboratory
% EPFL

clc; clear all; close all;

% Parameters
P.NumberOfFrames   = 5;
P.NumberOfBits     = 172;
P.Q_Ind = 12;

P.AccessType = 'CDMA';

P.ChannelType   = 'AWGN'; % 'Multipath', 'Fading', 'AWGN', 'ByPass'
if strcmp(P.ChannelType, 'Multipath') | strcmp(P.ChannelType, 'Fading')
    P.ChannelLength = 1; 
end
if strcmp(P.ChannelType, 'Fading')
    P.CoherenceTime = 19200/3; % A third of a second 
end

P.HadLen = 64; % Length of Hadamard Sequence

P.K = 9; % Length of convolutional encoder
P.ConvSeq = [753 561]; % Rate 1/2
P.Rate = length(P.ConvSeq);

P.LongCodeLength = 42; % PN Sequence

if strcmp(P.ChannelType, 'Multipath') | strcmp(P.ChannelType, 'Fading')
    P.RakeFingers = 1; 
    if (P.RakeFingers > P.ChannelLength)
        error('Fingers has to be smaller or equal to channels !')
    end
end

P.SequenceMask = [1,1,0,0,0,1,1,0,0,0, randi([0 1],1,32)];


P.SNRRange = -50:5:0; % SNR Range to simulate in dB

BER = SISOsimulator(P);

simlab = sprintf('%s - Length: %d - Users: %d' ,P.ChannelType,P.ChannelLength,P.CDMAUsers);

figure;
semilogy(P.SNRRange,BER,'b.-','DisplayName',simlab)

xlabel('SNR','FontSize',12,'FontWeight','bold');
ylabel('BER','FontSize',12,'FontWeight','bold');
xlim([min(P.SNRRange) max(P.SNRRange)]);
grid minor;
legend('-DynamicLegend');