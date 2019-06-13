% Wireless Receivers Project:
% Anael Buchegger, Tim Tuuva, David Sanchez
%
% SISO Simulation Framework
% 1 User only (simplest case) 
%
% Telecommunications Circuits Laboratory
% EPFL


function BER = SISOsimulator(P)
    
    % Generate the spreading sequence
    HadamardMatrix = hadamard(P.HadLen);     
    
    % Number of transmitted bits
    NbTXBits    = P.Rate*(P.NumberOfBits + P.Q_Ind + P.K-1);
    NumberOfChips  = NbTXBits*P.HadLen;
    
    HadamSequence     = HadamardMatrix(:,42); % 42 but could pick anything
    
    LongCode = comm.PNSequence('Polynomial',[42 7 6 5 3 2 1 0], ...
                               'Mask', P.SequenceMask, ...
                               'InitialConditions', randi([0 1],1,42), ...
                               'SamplesPerFrame', NbTXBits);

    % Generate a Code for each users                     
    PNSequence = step(LongCode).';
    
    % Channel
    switch P.ChannelType
        case {'Multipath', 'Fading'}
            NumberOfChipsRX = NumberOfChips+P.ChannelLength-1;
        otherwise
            NumberOfChipsRX = NumberOfChips;
    end

Results = zeros(1,length(P.SNRRange));

% Note that these are terminated, so they include the tail
trellis = poly2trellis(P.K, P.ConvSeq);
convEnc = comm.ConvolutionalEncoder(trellis, 'TerminationMethod', 'Terminated');
convDec = comm.ViterbiDecoder(trellis, 'TerminationMethod', 'Terminated', 'InputFormat','Hard');

for frame = 1:P.NumberOfFrames
    
    frame
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bits = randi([0 1],1,P.NumberOfBits); % Random Data
    
    % Add Frame Quality Indicator (bonus) Random as our data is random
    bits_ind = [bits randi([0 1],1,P.Q_Ind)];
    
    % Convolutional encoding
    encoded_bits = convEnc(bits_ind.').';
    
    % Interleaver
    if strcmp(P.Interleaving, 'On')
        encoded_bits=matintrlv(encoded_bits.',P.IntrlvRows,NbTXBits/P.IntrlvRows).';  
    end
    
    % Pulse Shape (PNSequence)
    PN_symbols = xor(PNSequence, encoded_bits);

    % Modulation : BPSK
    symbols = -(2*PN_symbols - 1);
    
    % Spreading with Hadamard
    symbol_spread = HadamSequence * symbols;    
    waveform = reshape(symbol_spread,1,P.HadLen*NbTXBits);
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Channel
    switch P.ChannelType
        case {'ByPass','AWGN'}
            himp = 1;
        case 'Multipath'
            himp = sqrt(1/2)* (randn(1,P.ChannelLength) + 1i * randn(1,P.ChannelLength));
            himp = himp/sqrt(sum(abs(himp).^2)); % Normalization
        case 'Fading'
            % Channel impulse for each path
            himp = sqrt(1/2)* (randn(1,P.ChannelLength) + 1i * randn(1,P.ChannelLength));
            himp = himp/sqrt(sum(abs(himp).^2)); % Normalization
            % Channel variation for each bit
            h = channel(P.ChannelLength,NumberOfChips,1,P.CoherenceTime,1);           
        otherwise
            disp('Channel not supported')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Simulation
    switch P.ChannelType
        case {'AWGN','ByPass'}
            snoise = randn(1,NumberOfChips) + 1i* randn(1,NumberOfChips) ;      
        case {'Multipath','Fading'}
            snoise = randn(1,NumberOfChipsRX) + ...
                      1i* randn(1,NumberOfChipsRX);           
    end
    
    % SNR Range
    for snr = 1:length(P.SNRRange)
        SNRdb  = P.SNRRange(snr);
        SNRlin = 10^(SNRdb/10);
        noise  = 1/sqrt(2*P.HadLen*SNRlin) *snoise;
        
        % Channel
        switch P.ChannelType
            case 'ByPass'
                y = waveform * himp;                
            case 'AWGN'
                y = waveform * himp + noise;              
            case 'Multipath'
                y = conv(waveform,himp) + noise;           
            case 'Fading'
                y = zeros(1,NumberOfChipsRX);
                y_conv = himp.'*waveform;
                for i = 1:P.ChannelLength
                    y_channel = y_conv(i,:) .* h(1,:,i);
                    y(i:i+NumberOfChips-1) = y(i:i+NumberOfChips-1) + y_channel;          
                end      
                y = y + noise;
            otherwise
                disp('Channel not supported')
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Receiver
        switch P.ChannelType
            case {'AWGN','ByPass'}
                % Despreading
                rxsymbols = conj(himp).*HadamSequence.'*reshape(y, P.HadLen, NbTXBits);
                desp_bits = rxsymbols < 0; 

            case {'Multipath','Fading'} % Rake receiver
                rxsymbols = zeros(P.RakeFingers, NbTXBits);
                
                % Order the best fingers
                [~,ind] = maxk(himp,P.RakeFingers);

                for finger = 1:P.RakeFingers
                    % Channel would be estimated for each path
                    himp_conj = conj(himp(ind(finger)));
                    
                    if strcmp(P.ChannelType,'Multipath')
                        % Despreading
                        rx_despread = HadamSequence.'*reshape(...
                                                      y(ind(finger):ind(finger)+NumberOfChips-1),...
                                                      P.HadLen, NbTXBits);
                    else % Fading
                        % Channel effect would be followed during transmission
                        h_conj = conj(h(1,:,ind(finger)));
                        % Despreading and neutralizing channel effect
                        rx_despread = HadamSequence.'*reshape(h_conj .*...
                                                      y(ind(finger):ind(finger)+NumberOfChips-1),...
                                                      P.HadLen, NbTXBits);
                    end
                    
                    % Neutralizing global channel effect
                    rxsymbols(finger,:) = himp_conj * rx_despread;
                end

                % Summing over all fingers to get some diversity
                desp_bits = sum(rxsymbols,1) < 0; 
                
            otherwise
                disp('Source Encoding not supported')
        end
        
        % UN-PN
        unPN_symbols = xor(PNSequence, desp_bits);
        
        % De-interleaver
        unPN_symbols=double(unPN_symbols);   
        if strcmp(P.Interleaving, 'On')
           unPN_symbols=matdeintrlv(unPN_symbols.',32,12).';
        end
        
        % Decoding Viterbi
        decoded_bits = convDec((unPN_symbols).').';
        
        % Remove the bit encoding trail
        rxbits = decoded_bits(:,1:end-P.Q_Ind-(P.K-1)); 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % BER count
        errors =  sum(sum(rxbits ~= bits));     
        Results(snr) = Results(snr) + errors;
        
    end
end

BER = Results/(P.NumberOfBits*P.NumberOfFrames);
end
