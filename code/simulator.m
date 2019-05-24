% Wireless Receivers Project:
% Anael Buchegger, Tim Tuuva, David Sanchez
%
% Direct Sequence Spread Spectrum Simulation Framework
%
% Telecommunications Circuits Laboratory
% EPFL


function BER = simulator(P)

    RX = P.RX_Users;
    
    % Generate the spreading sequences % Custom matrix here
    HadamardMatrix = hadamard(P.HadLen);    %TODO normalization      
    
    % 
    NumberOfBits   = P.NumberOfBits*P.Modulation; % per Frame
    
    %RATE_BITS = P.Rate * NumberOfBits;
    NbTXBits    = P.Rate*(P.NumberOfBits + P.Q_Ind + P.K-1);
    NumberOfChips  = NbTXBits*P.HadLen;
    
    HadamSequence     = HadamardMatrix(:,42);%genbarker(P.LongCodeLength); % -(2*step(GS)-1);
    
    LongCode = comm.PNSequence('Polynomial',[42 7 6 5 3 2 1 0], ...
                               'Mask', P.SequenceMask, ...
                               'InitialConditions', randi([0 1],1,42), ...
                               'SamplesPerFrame', NbTXBits);
    PNSequence = step(LongCode);
    
    % Channel
    switch P.ChannelType
        case 'Multipath'
            NumberOfChipsRX   = NumberOfChips+P.ChannelLength-1;
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

    bits = randi([0 1],1,NumberOfBits); % Random Data
    
    % Add Frame Quality Indicator (bonus)
    bits_ind = [bits randi([0 1],1,P.Q_Ind)];
    
    % Convolutional encoding
    encoded_bits = convEnc(bits_ind.');
    
    % Symbol repetition
    encoded_bits = repmat(encoded_bits, 1, NbTXBits/length(encoded_bits));
    encoded_bits = encoded_bits(:);
    
    % Here comes the interleaver (TODO)
    
    % Pulse Shape (PNSequence)    
    PN_symbols = xor(PNSequence, encoded_bits);

    % Modulation : BPSK
    symbols = -(2*PN_symbols - 1);
    
    % Spreading with Hadamard
    symbol_spread = HadamSequence * symbols.';
    waveform = reshape(symbol_spread, 1, P.HadLen*length(encoded_bits));    
 
    % Channel
    switch P.ChannelType
        case 'ByPass',
            himp = ones(P.ChannelLength,1,RX);
        case 'AWGN',
            himp = ones(P.ChannelLength,1,RX);
        case 'Multipath',
            himp = sqrt(1/2)* (randn(P.ChannelLength, RX) + 1i * randn(P.ChannelLength, RX));
        case 'Fading',
            himp = channel(P.ChannelLength,NumberOfChipsRX,1,P.CoherenceTime,1);
        
        otherwise,
            disp('Channel not supported')
    end
    
    mwaveform = repmat(waveform,[P.ChannelLength 1 RX]); %WTF happens here ?

    
    %%%
    % Simulation
    snoise = ( randn(1,NumberOfChips,RX) + 1i* randn(1,NumberOfChips,RX) );
    
    % SNR Range
    for ss = 1:length(P.SNRRange)
        SNRdb  = P.SNRRange(ss);
        SNRlin = 10^(SNRdb/10);
        noise  = 1/sqrt(2*P.HadLen*SNRlin) *snoise;
        
        % Channel
        switch P.ChannelType
            case 'ByPass',
                y = zeros(P.ChannelLength,NumberOfChipsRX,RX); %Normally add the users here!
                for r = 1:RX
                    for i = 1:P.ChannelLength
                        y(i,:,r) = conv(mwaveform(i,:,r),himp(i,:,r)); %TODO MIMO
                    end
                end
            case 'AWGN',
                y = zeros(P.ChannelLength,NumberOfChipsRX,RX); %Normally add the users here!
                for r = 1:RX
                    for i = 1:P.ChannelLength
                        y(i,:,r) = conv(mwaveform(i,:,r),himp(i,:,r)) + noise(i,:,r);  %TODO MIMO
                    end
                end
            case 'Multipath'     
                y = zeros(P.ChannelLength,NumberOfChipsRX,RX); %Normally add the users here!
                for i = 1:P.ChannelLength
                    y(i,i:NumberOfChips+i-1,RX) = conv(mwaveform(i,:,RX),himp(i,:)) + noise; 
                end
                
            case 'Fading',
                y = zeros(P.ChannelLength,NumberOfChipsRX,RX); %Normally add the users here!
                for i = 1:P.ChannelLength
                    y(i,:,RX) = mwaveform(i,:,RX) .* himp(RX,:,i) + noise;
                end
            otherwise,
                disp('Channel not supported')
        end
        
 
        % Receiver
        switch P.ReceiverType
            case 'Rake',  
                % Despreading
                rxsymbols = zeros(P.ChannelLength, NbTXBits, RX);
                
                for r = 1:RX
                    if ~strcmp(P.ChannelType,'Fading')
                        [~,ind] = maxk(himp,P.ChannelLength);
                    else
                        [himp_mean,ind] = maxk(mean(himp(RX,:,:)),P.ChannelLength);
                    end
                    for finger = 1:P.RakeFingers
                        if ~strcmp(P.ChannelType,'Fading')
                            
                            rxsymbols(finger,:,r) = conj(himp(ind(finger)))*HadamSequence.'*...
                                reshape(y(ind(finger),ind(finger):ind(finger)+NumberOfChips-1,RX), ...
                                P.HadLen, NumberOfChips/P.HadLen);
                        else
                            rxsymbols(finger,:,r) = conj(himp_mean(ind(finger)))*HadamSequence.'*...
                                reshape(y(ind(finger),:,RX), ...
                                P.HadLen, NumberOfChipsRX/P.HadLen);
                        end
                    end
                end
                % TODO is that ok??
                desp_bits = reshape(sum(rxsymbols,1) < 0,NbTXBits,RX); 
                
            otherwise,
                disp('Source Encoding not supported')
        end
        % UN-PN
        unpn_bits = xor(PNSequence, desp_bits);
        
        
        % Decoding Viterbi
        rxbits = convDec(double(unpn_bits));
        
        % Remove the 8 bit encoding trail
        rxbits = rxbits(1:end-P.Q_Ind-8).'; % TODO magick number

        % BER count
        errors =  sum(rxbits ~= bits);
        
        Results(ss) = Results(ss) + errors;
        
    end
end

BER = Results/(NumberOfBits*P.NumberOfFrames);
end