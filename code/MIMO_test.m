% Wireless Receivers Project:
% Anael Buchegger, Tim Tuuva, David Sanchez
%
% CDMA Parameter File
% MIMO configuration
%
% Telecommunications Circuits Laboratory
% EPFL

clc; clear all; close all;

% Parameters
P.NumberOfFrames      = 50;
P.NumberOfBits     = 172; 
P.Q_Ind = 12;

P.CDMAUsers     = 1; %TODO
P.RXperUser     = 2;
P.TXperUser     = 2; % TODO is that ok there isnot a equal number of RX and TX antennas?
                     % right now We can't have a different number, we have
                     % to decide which anteanns signals we just keep or mix
                     % in a good manner the information

P.ChannelType   = 'Multipath'; % 'Multipath', 'Fading', 'AWGN', 'ByPass'

if ((strcmp(P.ChannelType, 'ByPass') || strcmp(P.ChannelType, 'AWGN')) ...
        && P.RXperUser > 1 && P.TXperUser > 1)
    error('ByPass or AWGN make no sense with MIMO ! (destructive interferences)')
end

% Only applies for fading and multipath
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

% P.ReceiverType  = 'Rake'; % Only one type for the project

BER = MIMOsimulator(P);

if strcmp(P.ChannelType, 'Multipath') | strcmp(P.ChannelType, 'Fading')
    simlab = sprintf('%s - Paths: %d - TX/RX : %d/%d - Fingers : %d - Users: %d' ,...
         P.ChannelType,P.ChannelLength,P.TXperUser,P.RXperUser,P.RakeFingers,P.CDMAUsers);
else
    simlab = sprintf('%s - Users: %d' ,P.ChannelType,P.CDMAUsers);
end

figure;
semilogy(P.SNRRange,BER,'b.-','DisplayName',simlab)

xlabel('SNR','FontSize',12,'FontWeight','bold');
ylabel('BER','FontSize',12,'FontWeight','bold');
xlim([min(P.SNRRange) max(P.SNRRange)]);
grid minor;
legend('-DynamicLegend');