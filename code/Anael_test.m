% Wireless Receivers II - Assignment 2: TODO Rename to report?? :P
%
% Fading CDMA Parameter File
%
% Telecommunications Circuits Laboratory
% EPFL

clc; clear all; close all;

% Parameters
P.ReceiverType  = 'Rake';
P.NumberOfFrames      = 50;
P.NumberOfBits     = 172; % TODO, in simulator: NumberOfBits   = P.NumberOfSymbols*P.Modulation*RX; % per Frame
P.Q_Ind = 12;

P.AccessType = 'CDMA';
P.CDMAUsers     = 1;

P.Modulation    = 1;        % 1: BPSK, 2: 64ary

P.ChannelType   = 'Multipath'; % 'Multipath', 'Fading', 'AWGN', 'ByPass'
P.ChannelLength = 3;

P.HadLen = 64; % Length of Hadamard Sequence
P.HadIn = 6; % Length of Hadamard input (6 -> 64)

P.K = 9; % Length of convolutional encoder
P.ConvSeq = [753 561]; % Rate 1/2
P.Rate = length(P.ConvSeq);

P.LongCodeLength = 42; % PN Sequence
P.RakeFingers = 10; 

P.SequenceMask = [1,1,0,0,0,1,1,0,0,0, randi([0 1],1,32)];


P.SNRRange = -50:5:0; % SNR Range to simulate in dB

% P.ReceiverType  = 'Rake'; % Only one type for the project

BER = simulator(P);

simlab = sprintf('%s - Length: %d - Users: %d' ,P.ChannelType,P.ChannelLength,P.CDMAUsers);

figure(1)
semilogy(P.SNRRange,BER,'b.-','DisplayName',simlab)

xlabel('SNR','FontSize',12,'FontWeight','bold');
ylabel('BER','FontSize',12,'FontWeight','bold');
xlim([min(P.SNRRange) max(P.SNRRange)]);
grid minor;
legend('-DynamicLegend');