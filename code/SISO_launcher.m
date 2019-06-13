% Wireless Receivers Project:
% Anael Buchegger, Tim Tuuva, David Sanchez
%
% CDMA Parameter File
% SISO, 1 User only (simplest case)
%
% Telecommunications Circuits Laboratory
% EPFL

clc; clear all; %close all;

%% Parameters
P.NumberOfFrames = 10;
P.SNRRange = -50:1:-20; % SNR Range to simulate in dB

P.NumberOfBits = 172;
P.Q_Ind = 12;
P.HadLen = 64; % Length of Hadamard Sequence
P.K = 9; % Length of convolutional encoder
P.ConvSeq = [753 561]; % Rate 1/2
P.Rate = length(P.ConvSeq);
P.SequenceMask = [1,1,0,0,0,1,1,0,0,0, randi([0 1],1,32)]; %Mask for sequence
P.Interleaving = 'On'; % 'On' or 'Off'
P.IntrlvRows = 32; % For a rate of 9600bps

%% Channel and Detectors
P.ChannelType   = 'Multipath'; % 'ByPass','Multipath','Fading', 'AWGN' 

% Only applies for fading and multipath
P.ChannelLength = 3; 
P.RakeFingers = 3; 
% Check
if (P.RakeFingers > P.ChannelLength)
    error('Fingers has to be smaller or equal to channels !')
end

% Only applies for fading
P.CoherenceTime = 19; % A thousandth of a second 

%% Simulation
BER = SISOsimulator(P);

if strcmp(P.ChannelType, 'Multipath') | strcmp(P.ChannelType, 'Fading')
    simlab = sprintf('%s - Paths: %d - Fingers : %d' ,P.ChannelType,P.ChannelLength, P.RakeFingers);
else
    simlab = sprintf('%s' ,P.ChannelType);
end

%% Plotting
figure
semilogy(P.SNRRange,BER,'b.-','DisplayName',simlab)

xlabel('SNR','FontSize',12,'FontWeight','bold');
ylabel('BER','FontSize',12,'FontWeight','bold');
xlim([min(P.SNRRange) max(P.SNRRange)]);
grid minor;
legend('-DynamicLegend');