% Wireless Receivers Project:
% Anael Buchegger, Tim Tuuva, David Sanchez
%
% CDMA Parameter File
%
% Telecommunications Circuits Laboratory
% EPFL

clc; clear all; close all;

% Data parameters
P.ReceiverType  = 'Rake';
P.NumberOfFrames      = 5;
P.NumberOfBits     = 172;
P.Q_Ind = 12;

P.AccessType = 'CDMA';
P.CDMAUsers     = 1; %TODO implement this!
P.RXperUser     = 1;
P.TXperUser     = P.RXperUser; % TODO is that ok there isnot a equal number of RX and TX antennas?
                               % right now We can't have a different number, we have
                               % to decide which anteanns signals we just keep or mix
                               % in a good manner the information


P.Modulation    = 1;        % 1: BPSK, 2: 64ary

% Channel Modelisation
P.ChannelType   = 'AWGN'; % 'Multipath', 'Fading', 'AWGN', 'ByPass'
P.ChannelLength = 1;  %TODO bug, Channellength cant be bigger than RX ;-/ is that normal??
P.RakeFingers = 1; 
P.CoherenceTime = 19200/3; % A third of a second

if (P.RakeFingers > P.ChannelLength)
    error('Fingers has to be smaller or equal to channels !')
end

P.SNRRange = -50:5:0; % SNR Range to simulate in dB

% Spreading Walsh function
P.HadLen = 64; % Length of Hadamard Sequence

% Conv encodr
P.K = 9; % Length of convolutional encoder
P.ConvSeq = [753 561]; % Rate 1/2
P.Rate = length(P.ConvSeq);

% PN Sequence
P.LongCodeLength = 42; 
P.SequenceMask = [1,1,0,0,0,1,1,0,0,0, randi([0 1],1,32)];

% Simulation
BER = simulator(P);

simlab = sprintf('%s - Length: %d - Users: %d' ,P.ChannelType,P.ChannelLength,P.CDMAUsers);

figure;
semilogy(EbN0shift(P),BER,'b.-','DisplayName',simlab)

xlabel('EbN0','FontSize',12,'FontWeight','bold');
ylabel('BER','FontSize',12,'FontWeight','bold');
xlim([min(EbN0shift(P)) max(EbN0shift(P))]);
grid minor;
legend('-DynamicLegend');