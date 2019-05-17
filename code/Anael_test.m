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
P.NumberOfSymbols     = 172; % TODO, in simulator: NumberOfBits   = P.NumberOfSymbols*P.Modulation*RX; % per Frame

P.AccessType = 'CDMA';
P.CDMAUsers     = 1;

P.Modulation    = 1;        % 1: BPSK, 2: 64ary

P.ChannelType   = 'AWGN'; % 'Multipath', 'Fading'
P.ChannelLength = 1;

P.HamLen = 64; % Length of Hadamard Sequence
P.HadIn = 6; % Length of Hadamard input (6 -> 64)

P.K = 9; % Length of convolutional encoder
P.ConvSeq = [557 663 711]; % Rate 1/3
P.Rate = 3;

P.LongCodeLength = 42; % PN Sequence
P.RakeFingers = 2; 

P.SNRRange = -50:10:0; % SNR Range to simulate in dB

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