% Wireless Receivers Project:
% Anael Buchegger, Tim Tuuva, David Sanchez
%
% SISO, 1 User only (simplest case)
%
% Telecommunications Circuits Laboratory
% EPFL


function BER = SISOsimulator(P)
    
    % Generate the spreading sequence
    HadamardMatrix = hadamard(P.HadLen)/sqrt(P.HadLen);    %TODO normalization      
    
    NbTXBits    = P.Rate*(P.NumberOfBits + P.Q_Ind + P.K-1);
    NumberOfChips  = NbTXBits*P.HadLen;
    
    HadamSequence     = HadamardMatrix(:,42); %42 but could pick anything
    
    LongCode = comm.PNSequence('Polynomial',[42 7 6 5 3 2 1 0], ...
                               'Mask', P.SequenceMask, ...
                               'InitialConditions', randi([0 1],1,42), ...
                               'SamplesPerFrame', NbTXBits);
                           
    PNSequence = step(LongCode).';
    
    % Channel
    switch P.ChannelType
        case {'Multipath', 'Fading'}
            NumberOfChipsRX = NumberOfChips+P.ChannelLength-1;
        otherwise
            NumberOfChipsRX = NumberOfChips;
    end

Results = zeros(1,length(P.SNRRange));

% verifier si c'est les bon coeff!
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
    encoded_bits = convEnc(bits_ind.').'; % TODO doesnt give the same encoding signal, why??? IMPORTANT
    
    % Here comes the interleaver (TODO)
    
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
        case 'Multipath',
            himp = sqrt(1/2)* (randn(1,P.ChannelLength) + 1i * randn(1,P.ChannelLength));
        case 'Fading',
            himp = channel(P.ChannelLength,NumberOfChips,1,P.CoherenceTime,1);           
        otherwise,
            disp('Channel not supported')
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Simulation
    
    switch P.ChannelType
        case {'AWGN','ByPass'}
            snoise = randn(1,NumberOfChips) + 1i* randn(1,NumberOfChips) ;
            
        case {'Multipath','Fading'}
            snoise = randn(P.ChannelLength,NumberOfChips) + ...
                      1i* randn(P.ChannelLength,NumberOfChips);    
            mwaveform = repmat(waveform,[P.ChannelLength 1]);
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
                y = zeros(P.ChannelLength,NumberOfChipsRX); %Normally add the users here!
                for i = 1:P.ChannelLength
                    y(i,i:i+NumberOfChips-1) = conv(mwaveform(i,:),himp(i)) + noise(i,:);
                end               
            case 'Fading'
                y = zeros(P.ChannelLength,NumberOfChipsRX);
                for i = 1:P.ChannelLength
                    y(i,i:i+NumberOfChips-1) = mwaveform(i,:) .* himp(1,:,i) + noise(i,:);
                end
            otherwise
                disp('Channel not supported')
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Receiver
        switch P.ChannelType % Rake receiver makes no sense in this case
            case {'AWGN','ByPass'}
                % Despreading
                rxsymbols = conj(himp).*HadamSequence.'*reshape(y, P.HadLen, NbTXBits);
                    
                desp_bits = rxsymbols < 0; 

            case {'Multipath','Fading'} % Rake receiver
                rxsymbols = zeros(P.RakeFingers, NbTXBits);
                
                % Order the best fingers
                if strcmp(P.ChannelType,'Multipath')
                    [~,ind] = maxk(himp,P.RakeFingers);
                else
                    % Select the best fingers based on mean channel
                    % (estimated ?)
                    [~,ind] = maxk(mean(himp(1,:,:),2),P.RakeFingers);
                end        

                for finger = 1:P.RakeFingers
                    if strcmp(P.ChannelType,'Multipath')
                        % Channel would be estimated for each path
                        rx_channel = conj(himp(ind(finger)));
                        % Despreading
                        rx_despread = HadamSequence.'*reshape(...
                                                      y(ind(finger),ind(finger):ind(finger)+NumberOfChips-1),...
                                                      P.HadLen, NbTXBits);
                        % Neutralizing channel effect
                        rxsymbols(finger,:) = rx_channel * rx_despread;
                    else
                        % Channel would be estimated for each path, and
                        % followed during transmission
                        rx_channel = conj(himp(1,:,ind(finger)));
                        % Despreading and neutralizing channel effect
                        rxsymbols(finger,:) = HadamSequence.'*reshape(rx_channel .*...
                                                      y(ind(finger),ind(finger):ind(finger)+NumberOfChips-1),...
                                                      P.HadLen, NbTXBits);
                    end
                end

                % Summing over all fingers to get some diversity
                desp_bits = sum(reshape(sum(rxsymbols,1) < 0,1,NbTXBits),1); 
                
            otherwise,
                disp('Source Encoding not supported')
        end
        % UN-PN
        unPN_symbols = xor(PNSequence, desp_bits);
        
        % Decoding Viterbi
        decoded_bits = convDec(double(unPN_symbols).').'; % TODO, beurk beurk no??
        
        % Remove the bit encoding trail
        rxbits = decoded_bits(:,1:end-P.Q_Ind-(P.K-1)); 
        
        % Here comes the de-interleaver (TODO)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % BER count
        errors =  sum(sum(rxbits ~= bits));
        
        Results(snr) = Results(snr) + errors;
        
    end
end

BER = Results/(P.NumberOfBits*P.NumberOfFrames);
end