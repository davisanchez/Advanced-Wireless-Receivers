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
    NumberOfChips  = 24576;%P.NumberOfBits *P.Modulation*P.HadLen/P.HadIn*P.LongCodeLength*P.Rate;
    % NumberOfChips  = P.HadLen * RATE_BITS/P.HadIn; % per Frame

    HadamSequence     = HadamardMatrix(:,42);%genbarker(P.LongCodeLength); % -(2*step(GS)-1);
    
    LongCode = comm.PNSequence('Polynomial',[42 7 6 5 3 2 1 0], ...
                               'Mask', P.SequenceMask, ...
                               'InitialConditions', randi([0 1],1,42), ...
                               'SamplesPerFrame', 384);
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
convDec = comm.ViterbiDecoder(trellis, 'TerminationMethod', 'Terminated');

for ii = 1:P.NumberOfFrames
    
    ii

    bits = randi([0 1],1,NumberOfBits); % Random Data
    
    % TODO Add Frame Quality Indicator (bonus)
    bits = [bits randi([0 1],1,P.Q_Ind)];
    
    % TODO: Convolutional encoding
    encoded_bits = step(convEnc,bits.');
    %encoded_bits = convenc(bits,trellis);
    % TODO taking only the first stream?? dunno we got 540 and 540/3=180
    % which is 172+8 bits! No we take everything, we increase the rate
    % conv_bits = conv_bits(1:NumberOfBits);
    
    % Symbol repetition
    encoded_bits = repmat(encoded_bits, 1, 384/length(encoded_bits)); %TODO 384 faire une variable avec ca
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
            himp = sqrt(1/2)* ( randn(RX,P.ChannelLength) + 1i * randn(RX,P.ChannelLength) );
%             himp = (ones(RX,1) * sqrt(P.PDP)) .* himp;
        otherwise,
            disp('Channel not supported')
    end
    
    mwaveform = repmat(waveform,[P.RakeFingers 1 RX]); %WTF happens here ?

    
    %%%
    % Simulation
    snoise = ( randn(1,length(waveform),RX) + 1i* randn(1,length(waveform),RX) );
    
    % SNR Range
    for ss = 1:length(P.SNRRange)
        SNRdb  = P.SNRRange(ss);
        SNRlin = 10^(SNRdb/10);
        noise  = 1/sqrt(2*P.HadLen*SNRlin) *snoise;
        
        % Channel
        switch P.ChannelType
            case 'ByPass',
                y = zeros(P.RakeFingers,length(waveform),RX); %Normally add the users here!
                for i = 1:P.RakeFingers
                    y(i,:,RX) = conv(mwaveform(i,:,RX),himp(i,:)); 
                end
            case 'AWGN',
                y = zeros(P.RakeFingers,length(waveform),RX); %Normally add the users here!
                for i = 1:P.RakeFingers
                    y(i,:,RX) = conv(mwaveform(i,:,RX),himp(i,:)) + noise; 
                end
            case 'Multipath'     
                y = zeros(1,NumberOfChips+P.ChannelLength-1,RX);
                for i = 1:RX
                    y(1,:,i) = conv(mwaveform(1,:,i),himp(i,:)) + noise(1,:,i); 
                end
            otherwise,
                disp('Channel not supported')
        end
        
 
        % Receiver
        switch P.ReceiverType
            case 'Rake',  
                % Despreading
                [~,ind] = maxk(himp,P.RakeFingers);
                for(finger = 1:P.RakeFingers)
                    rxsymbols(finger,:) = conj(himp(ind(finger)))*HadamSequence.'*...
                                          reshape(y(ind(finger),1:NumberOfChips,RX), ...
                                          P.HadLen, NumberOfChips/P.HadLen);
                end
                desp_bits = (0 > real(rxsymbols(1,:))).'; %TODO correct this%reshape(sum(rxsymbols,1) < 0,1,P.NumberOfSymbols);
                
%                 % Hadamard
%                 for j=1:RX
%                     y_rx=rxbits;%sign(real(y(:,:,j))); %TODO hard decision, good?
% 
%                     for i=1:P.ChannelLength    
%                         %TODO reshape(y_rx(i:i+NumberOfChips-1),SeqLen,NumberOfBits/RX); 
%                         
%                         
%                         % TODO remove this???
%                         %rxsymbols(i,:)=SpreadSequence(:,j).'*y_reshape;
%                         %rxsymbols(i,:)=rxsymbols(i,:)*conj(himp(j,i));
%                     end
%                     %TODO: reshape(sum(rxsymbols,1) < 0,1,P.NumberOfSymbols);
%                     %rxbits(j:RX:RX*RATE_BITS) = reshape(sum(rxsymbols,1) < 0,1,RATE_BITS);
%                 end
            otherwise,
                disp('Source Encoding not supported')
        end
        % UN-PN
        unpn_bits = xor(PNSequence, desp_bits);
        
        
        % TODO: Decoding Viterbi, is that correct??? Dunno if we have to use it
        % :S
        %rxbits = step(convDec,unpn_bits); % right now it doesnt work,
        % problem of shape and size of array
        rxbits = vitdec(unpn_bits,trellis,34,'trunc','hard');
        
        % Remove the 8 bit encoding trail ????? TODO
        rxbits = rxbits(1:end-8).'; % TODO magick number

        % BER count
        errors =  sum(rxbits ~= bits);
        
        Results(ss) = Results(ss) + errors;
        
    end
end

BER = Results/(NumberOfBits*P.NumberOfFrames);
end

function seq = genbarker(len)
    BarkerSeq = [+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1];

    factor = ceil(len/length(BarkerSeq));
    b = repmat(BarkerSeq,1,factor);
    b = BarkerSeq.'*ones(1,factor);
    seq = b(1:len).';
end
