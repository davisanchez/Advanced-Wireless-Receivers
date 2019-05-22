% Wireless Receivers Project:
% Anael Buchegger, Tim Tuuva, David Sanchez
%
% Direct Sequence Spread Spectrum Simulation Framework
%
% Telecommunications Circuits Laboratory
% EPFL


function BER = simulator(P)

    RX = P.CDMAUsers;
    
    % Generate the spreading sequences % Custom matrix here
    HadamardMatrix = hadamard(P.HadLen);%/sqrt(P.HadLen);   TODO          
    
    %TODO magick number 3x172 -> 1/rate X Number of bits    
    NumberOfBits   = P.NumberOfBits*P.Modulation*RX; % per Frame
    
    RATE_BITS = P.Rate * NumberOfBits;
    NbTXBits    = P.Rate*(P.NumberOfBits + P.Q_Ind + 8);
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

for ii = 1:P.NumberOfFrames
    
    ii

    bits = randi([0 1],1,NumberOfBits); % Random Data
    
    % TODO Add Frame Quality Indicator (bonus)
    bits_ind = [bits randi([0 1],1,P.Q_Ind)];
    
    % TODO: Convolutional encoding
    encoded_bits = convEnc(bits_ind.');
    
    % Symbol repetition
    encoded_bits = repmat(encoded_bits, 1, NbTXBits/length(encoded_bits));
    encoded_bits = encoded_bits(:);
    
    % Here comes the interleaver (TODO)
    
    % Pulse Shape (PNSequence)    
    PN_symbols = xor(PNSequence, encoded_bits);

    % Modulation : Modulated with the 64ary
    switch P.Modulation % Modulate Symbols
        case 1 % BPSK
            symbols = -(2*PN_symbols - 1);
        case 2 % 64ary
%             symbols = reshape(conv_bits,P.HadIn,length(conv_bits)/P.HadIn);
%             had_index = bi2de(symbols.'); % Not sure...
%             symbols = HadamardMatrix(had_index(:)+1,:);
        otherwise
            disp('Modulation not supported')
    end
    
    
    % TODO: Spreading with Hadamard
    symbol_spread = HadamSequence * symbols.';
    waveform = reshape(symbol_spread, 1, P.HadLen*length(encoded_bits));    

    
    % Channel
    switch P.ChannelType
        case 'ByPass',
            himp = ones(P.RakeFingers,1);
        case 'AWGN',
            himp = ones(P.RakeFingers,1);
        case 'Multipath',
            himp = sqrt(1/2)* (randn(P.RakeFingers,P.ChannelLength) + 1i * randn(P.RakeFingers,P.ChannelLength));
        otherwise,
            disp('Channel not supported')
    end
    
    mwaveform = repmat(waveform,[P.RakeFingers 1 RX]); %WTF happens here ?

    
    %%%
    % Simulation
    snoise = ( randn(1,NumberOfChipsRX,RX) + 1i* randn(1,NumberOfChipsRX,RX) );
    
    % SNR Range
    for ss = 1:length(P.SNRRange)
        SNRdb  = P.SNRRange(ss);
        SNRlin = 10^(SNRdb/10);
        noise  = 1/sqrt(2*P.HadLen*SNRlin) *snoise;
        
        % Channel
        switch P.ChannelType
            case 'ByPass',
                y = zeros(P.RakeFingers,NumberOfChipsRX,RX); %Normally add the users here!
                for i = 1:P.RakeFingers
                    y(i,:,RX) = conv(mwaveform(i,:,RX),himp(i,:)); 
                end
            case 'AWGN',
                y = zeros(P.RakeFingers,NumberOfChipsRX,RX); %Normally add the users here!
                for i = 1:P.RakeFingers
                    y(i,:,RX) = conv(mwaveform(i,:,RX),himp(i,:)) + noise; 
                end
            case 'Multipath'     
                y = zeros(P.RakeFingers,NumberOfChipsRX+P.RakeFingers,RX); %Normally add the users here!
                for i = 1:P.RakeFingers
                    y(i,i:NumberOfChipsRX+i-1,RX) = conv(mwaveform(i,:,RX),himp(i,:)) + noise; 
                end
            otherwise,
                disp('Channel not supported')
        end
        
 
        % Receiver
        switch P.ReceiverType
            case 'Rake',  
                % Despreading
                rxsymbols = zeros(P.RakeFingers, NbTXBits);
                [~,ind] = maxk(himp,P.RakeFingers);
                for finger = 1:P.RakeFingers
                    rxsymbols(finger,:) = conj(himp(ind(finger)))*HadamSequence.'*...
                                          reshape(y(ind(finger),ind(finger):ind(finger)+NumberOfChips-1,RX), ...
                                          P.HadLen, NumberOfChips/P.HadLen);
                end
                desp_bits = reshape(sum(rxsymbols,1) < 0,1,NbTXBits).'; 
                
            otherwise,
                disp('Source Encoding not supported')
        end
        % UN-PN
        unpn_bits = xor(PNSequence, desp_bits);
        
        
        % TODO: Decoding Viterbi
        rxbits = convDec(double(unpn_bits));
        
        % Remove the 8 bit encoding trail ????? TODO
        rxbits = rxbits(1:end-P.Q_Ind-8).'; % TODO magick number

        % BER count
        errors =  sum(rxbits ~= bits);
        
        Results(ss) = Results(ss) + errors;
        
    end
end

BER = Results/(NumberOfBits*P.NumberOfFrames);
end