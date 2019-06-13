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
P.NumberOfFrames = 10;
P.SNRRange = -50:10:20;         % SNR Range to simulate in dB
P.NumberOfBits = 172;           % Number of transmited bits
P.Q_Ind = 12;                   % Extra bits for encoding
P.HadLen = 64;                  % Length of Hadamard Sequence
P.K = 9;                        % Length of convolutional encoder
P.ConvSeq = [753 561];          % Convolutionnal polynome with rate 1/2
P.Rate = length(P.ConvSeq);
                                %Mask for sequence
P.SequenceMask = [1,1,0,0,0,1,1,0,0,0, randi([0 1],1,32)];
P.Interleaving = 'On';          % 'On' or 'Off'
P.IntrlvRows = 32;              % For a rate of 9600bps
P.Decision='Soft';              % Hard or Soft decision for Viterbi decoder
%% Users and Antennas
P.CDMAUsers     = 4;            % Number of CDMA users
P.RXperUser     = 2;            % Rx antenna for each Users
P.TXperUser     = 2;            % Tx antenna for each Users

%% Channel and Detectors
P.ChannelLength = 3;           
P.RakeFingers = 2; 
P.Detector = 'ZF';              % 'ZF', 'MMSE', 'SIC'
P.Mode = 'HighDiversity';       % 'HighRate' or 'HighDiversity'

%% Checks
if(strcmp(P.Detector, 'ZF') && P.TXperUser > P.RXperUser)
    error('Can not have nTx > nRx, with ZF detector');
end

if (P.RakeFingers > P.ChannelLength)
    error('Fingers has to be smaller or equal to channels !')
end

if (strcmp(P.Detector, 'SIC') && strcmp(P.Mode,'HighDiversity'))
    error('High diversity mode has not been implemented with SIC !')
end    

%% Simulation
BER = MIMOsimulator(P);

simlab = sprintf('%s, %s Detector - TX/RX : %d/%d - Users: %d\nPaths: %d - Fingers : %d - %s Decision' ,...
         'Multipath',P.Detector,P.TXperUser,P.RXperUser,P.CDMAUsers,P.ChannelLength,P.RakeFingers,P.Decision);

%% Plotting
figure
semilogy(P.SNRRange,BER,'b.-','DisplayName',simlab)
title([P.Mode,' mode, ',num2str(P.NumberOfFrames),' frames, interleaver ',P.Interleaving]);
xlabel('SNR','FontSize',12,'FontWeight','bold');
ylabel('BER','FontSize',12,'FontWeight','bold');
xlim([min(P.SNRRange) max(P.SNRRange)]);
grid minor;
legend('-DynamicLegend');