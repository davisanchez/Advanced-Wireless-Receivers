% Wireless Receivers Project:
% Anael Buchegger, Tim Tuuva, David Sanchez
%
% CDMA Parameter File
%
% Telecommunications Circuits Laboratory
% EPFL

clc; clear all; close all;

% Parameters
%P.ReceiverType  = 'Rake';
P.NumberOfFrames      = 5;
P.NumberOfBits     = 172; % TODO, in simulator: NumberOfBits   = P.NumberOfSymbols*P.Modulation*RX; % per Frame
P.Q_Ind = 12;

P.AccessType = 'CDMA';
P.CDMAUsers     = 1; %TODO
P.RXperUser     = 2;
P.TXperUser     = 2; % TODO is that ok there isnot a equal number of RX and TX antennas?
                     % right now We can't have a different number, we have
                     % to decide which anteanns signals we just keep or mix
                     % in a good manner the information


P.Modulation    = 1;        % 1: BPSK, 2: 64ary

P.ChannelType   = 'ByPass'; % 'Multipath', 'Fading', 'AWGN', 'ByPass'
P.ChannelLength = 1;
P.CoherenceTime = 19200/3; % A third of a second

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

% P.ReceiverType  = 'Rake'; % Only one type for the project

BER = simulator(P);

simlab = sprintf('%s - Length: %d - Users: %d' ,P.ChannelType,P.ChannelLength,P.CDMAUsers);

figure;
semilogy(P.SNRRange,BER,'b.-','DisplayName',simlab)

xlabel('SNR','FontSize',12,'FontWeight','bold');
ylabel('BER','FontSize',12,'FontWeight','bold');
xlim([min(P.SNRRange) max(P.SNRRange)]);
grid minor;
legend('-DynamicLegend');