% Wireless Receivers Project:
% Anael Buchegger, Tim Tuuva, David Sanchez
%
% CDMA Parameter File
% MIMO configuration
%
% Telecommunications Circuits Laboratory
% EPFL

clc; clear all; % close all;

%% Parameters
P.NumberOfFrames = 100;
P.SNRRange = -50:1:20; % SNR Range to simulate in dB

P.NumberOfBits = 172; 
P.Q_Ind = 12;
P.HadLen = 64; % Length of Hadamard Sequence
P.K = 9; % Length of convolutional encoder
P.ConvSeq = [753 561]; % Rate 1/2
P.Rate = length(P.ConvSeq);
P.SequenceMask = [1,1,0,0,0,1,1,0,0,0, randi([0 1],1,32)]; %Mask for sequence
P.Interleaving = 'On'; % 'On' or 'Off'

%% Users and Antennas
P.CDMAUsers     = 1; %TODO
P.RXperUser     = 2;
P.TXperUser     = 2;

%% Channel and Detectors
P.ChannelType   = 'Multipath'; % 'Multipath', 'AWGN', 'ByPass'

% Only applies for multipath
P.ChannelLength = 2; 
P.RakeFingers = 1; 
P.Detector = 'ZF'; % 'ZF', 'MMSE', 'SIC'
P.Mode = 'HighRate'; % 'HighRate' or 'HighDiversity'

%% Checks
if(strcmp(P.Detector, 'ZF') && P.TXperUser > P.RXperUser)
    error('Can not have nTx > nRx, with ZF detector');
end

if ((strcmp(P.ChannelType, 'ByPass') || strcmp(P.ChannelType, 'AWGN')) ...
        && (P.RXperUser > 1 || P.TXperUser > 1))
    error('ByPass or AWGN make no sense with MIMO ! (destructive interferences)')
end

if (P.RakeFingers > P.ChannelLength)
    error('Fingers has to be smaller or equal to channels !')
end

if (strcmp(P.Detector, 'SIC') && strcmp(P.Mode,'HighDiversity'))
    error('High diversity mode has not been implemented with SIC !')
end    

%% Simulation
BER = MIMOsimulator(P);

if strcmp(P.ChannelType, 'Multipath')
    simlab = sprintf('%s, %s Detector - TX/RX : %d/%d\nPaths: %d - Fingers : %d \n%s mode - Users: %d' ,...
         P.ChannelType,P.Detector,P.TXperUser,P.RXperUser,P.ChannelLength,P.RakeFingers,P.Mode,P.CDMAUsers);
else
    simlab = sprintf('%s - Users: %d' ,P.ChannelType,P.CDMAUsers);
end

%% Plotting
figure
semilogy(P.SNRRange,BER,'b.-','DisplayName',simlab)

xlabel('SNR','FontSize',12,'FontWeight','bold');
ylabel('BER','FontSize',12,'FontWeight','bold');
xlim([min(P.SNRRange) max(P.SNRRange)]);
grid minor;
legend('-DynamicLegend');