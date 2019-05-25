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
    % TODO MIMO ??? is that ok?
    PNSequence = zeros(RX, NbTXBits);
    for r = 1:RX
        PNSequence(r, :) = step(LongCode);
    end
    
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

    bits = randi([0 1],RX,NumberOfBits); % Random Data
    
    % TODO debug MIMO
    %bits = randi([0 1],1,NumberOfBits); % Random Data
    %bits = [bits; bits];
    
    % Add Frame Quality Indicator (bonus)
    bits_ind = [bits randi([0 1],RX,P.Q_Ind)];
    
    % Convolutional encoding
    encoded_bits = zeros(RX, NbTXBits);
    for r=1:RX
        encoded_bits(r,:) = convEnc(bits_ind(r,:).'); % TODO doesnt give the same encoding signal, why??? IMPORTANT
    end
    
    % Symbol repetition
    encoded_bits = repmat(encoded_bits, 1, NbTXBits/length(encoded_bits));
    %encoded_bits = encoded_bits(:); TODO useless??
    
    % Here comes the interleaver (TODO)
    
    % Pulse Shape (PNSequence)    
    PN_symbols = xor(PNSequence, encoded_bits);

    % Modulation : BPSK
    symbols = -(2*PN_symbols - 1);
    
    % Spreading with Hadamard
    symbol_spread = zeros(RX, P.HadLen, NbTXBits);
    for r = 1:RX
        symbol_spread(r,:,:) = HadamSequence * symbols(r,:); % TODO MIMO, what shall we do? Different Hadamard sequence?
    end
    waveform = reshape(symbol_spread, RX, P.HadLen*NbTXBits);    
 
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
                        [himp_mean,ind] = maxk(mean(himp(r,:,:)),P.ChannelLength);
                    end
                    for finger = 1:P.RakeFingers
                        if ~strcmp(P.ChannelType,'Fading')
                            
                            rxsymbols(finger,:,r) = conj(himp(ind(finger)))*HadamSequence.'*...
                                reshape(y(ind(finger),ind(finger):ind(finger)+NumberOfChips-1,r), ...
                                P.HadLen, NumberOfChips/P.HadLen);
                        else
                            rxsymbols(finger,:,r) = conj(himp_mean(ind(finger)))*HadamSequence.'*...
                                reshape(y(ind(finger),:,r), ...
                                P.HadLen, NumberOfChipsRX/P.HadLen);
                        end
                    end
                end
                % TODO diversity or High rate possibility???
                desp_bits = reshape(sum(rxsymbols,1) < 0,NbTXBits,RX).'; 
                
            otherwise,
                disp('Source Encoding not supported')
        end
        % UN-PN
        unpn_bits = xor(PNSequence, desp_bits);
        
        
        % Decoding Viterbi
        decoded_bits = zeros(RX, NbTXBits/2);
        for r = 1:RX
            decoded_bits(r,:) = convDec(double(unpn_bits(r,:)).').'; % TODO, beurk beurk no??
        end
        
        % Remove the 8 bit encoding trail
        rxbits = decoded_bits(:,1:end-P.Q_Ind-8); % TODO magick number

        % BER count
        errors =  sum(sum(rxbits ~= bits)); % TODO good way to compute error here?
        
        Results(ss) = Results(ss) + errors;
        
    end
end

BER = Results/(NumberOfBits*P.NumberOfFrames);
end