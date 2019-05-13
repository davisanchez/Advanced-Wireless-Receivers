% Wireless Receivers II - Assignment 2: TODO Rename to report?? :P
%
% Fading CDMA Parameter File
%
% Telecommunications Circuits Laboratory
% EPFL

% Parameters
P.NumberOfFrames      = 500;
P.NumberOfSymbols     = 172; % TODO, in simulator: NumberOfBits   = P.NumberOfSymbols*P.Modulation*RX; % per Frame

P.AccessType = 'CDMA';
P.CDMAUsers     = 1;

P.Modulation    = 1;        % 1: BPSK

P.ChannelType   = 'AWGN'; % 'Multipath', 'Fading'
P.ChannelLength = 1;

P.HamLen = 8; % Length of Hadamard Sequence

P.SNRRange = -10:20; % SNR Range to simulate in dB

P.ReceiverType  = 'Rake';

BER = simulator(P);

simlab = sprintf('%s - Length: %d - Users: %d' ,P.ChannelType,P.ChannelLength,P.CDMAUsers);

figure(1)
semilogy(P.SNRRange,BER,'b.-','DisplayName',simlab)

xlabel('SNR','FontSize',12,'FontWeight','bold');
ylabel('BER','FontSize',12,'FontWeight','bold');
xlim([min(P.SNRRange) max(P.SNRRange)]);
grid minor;
legend('-DynamicLegend');