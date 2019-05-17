% Wireless Receivers II - Assignment 2:
%
% Direct Sequence Spread Spectrum Simulation Framework
%
% Telecommunications Circuits Laboratory
% EPFL


function BER = simulator(P)

    RX = P.CDMAUsers;
    
    % Generate the spreading sequences % Custom matrix here
    HadamardMatrix = hadamard(P.HamLen);%/sqrt(P.HamLen);   TODO          
    SpreadSequence = HadamardMatrix;
    
    SeqLen         = P.HamLen; 
    
    %TODO magick number 3x172 -> 1/rate X Number of bits    
    NumberOfBits   = P.NumberOfSymbols*P.Modulation*RX; % per Frame
    
    RATE_BITS = P.Rate * NumberOfBits;
    NumberOfChips  = P.NumberOfSymbols*P.Modulation*SeqLen*P.LongCodeLength*P.Rate/P.HadIn;
    % NumberOfChips  = P.HamLen * RATE_BITS/P.HadIn; % per Frame

    PNSequence     = genbarker(P.LongCodeLength); % -(2*step(GS)-1);
    
    
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
    
    % TODO: Convolutional encoding
    % conv_bits = step(convEnc,bits.');
    encoded_bits = convenc(bits,trellis);
    % TODO taking only the first stream?? dunno we got 540 and 540/3=180
    % which is 172+8 bits! No we take everything, we increase the rate
    % conv_bits = conv_bits(1:NumberOfBits);
 
    % Here comes the interleaver (TODO)
    
    
    % Modulation : Modulated with the 64ary
    switch P.Modulation % Modulate Symbols
        case 1 % BPSK
            symbols = -(2*encoded_bits - 1);
        case 2 % 64ary
%             symbols = reshape(conv_bits,P.HadIn,length(conv_bits)/P.HadIn);
%             had_index = bi2de(symbols.'); % Not sure...
%             symbols = HadamardMatrix(had_index(:)+1,:);
        otherwise
            disp('Modulation not supported')
    end
    

    % distribute symbols on users
    %TODO length(symbols) <- NumberOfBits
    SymUsers = reshape(symbols,RX,length(symbols)/RX);
    
    % TODO: Upsampling with Hadamard
    hada_bits = HadamardMatrix(bi2de(reshape(encoded_bits, length(encoded_bits)/P.HadIn, P.HadIn))+1,:);
    waveform = reshape(hada_bits, 1, P.HamLen*length(encoded_bits)/P.HadIn);
    
    % Pulse Shape (PNSequence)
    spread_waveform = PNSequence*waveform;
    waveform  = reshape(spread_waveform,1,NumberOfChips);

    mwaveform = repmat(waveform,[P.RakeFingers 1 RX]); %WTF happens here ?
    
    % Channel
    switch P.ChannelType
        case 'AWGN',
            himp = ones(P.RakeFingers,1);
        case 'Multipath',
            himp = sqrt(1/2)* ( randn(RX,P.ChannelLength) + 1i * randn(RX,P.ChannelLength) );
%             himp = (ones(RX,1) * sqrt(P.PDP)) .* himp;
        otherwise,
            disp('Channel not supported')
    end
    
    %%%
    % Simulation
    snoise = ( randn(1,NumberOfChipsRX,RX) + 1i* randn(1,NumberOfChipsRX,RX) );
    
    % SNR Range
    for ss = 1:length(P.SNRRange)
        SNRdb  = P.SNRRange(ss);
        SNRlin = 10^(SNRdb/10);
        noise  = 1/sqrt(2*SeqLen*SNRlin) *snoise;
        
        % Channel
        switch P.ChannelType
            case 'AWGN',
                y = zeros(P.RakeFingers,NumberOfChips,RX); %Normally add the users here!
                for i = 1:P.RakeFingers
                    y(i,:,RX) = conv(mwaveform(i,:,RX),himp(i,:));% + noise(1,:,i); 
                end
                %y = mwaveform;% + noise;
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
                    rxsymbols(finger,:) = conj(himp(ind(finger)))*PNSequence.'*...
                        reshape(y(ind(finger),1:NumberOfChips,RX),...
                                P.LongCodeLength, NumberOfChips/P.LongCodeLength)/42;
                end
                rxbits = rxsymbols(1,:); %TODO correct this%reshape(sum(rxsymbols,1) < 0,1,P.NumberOfSymbols);
                
                % Hadamard
                for j=1:RX
                    y_rx=sign(real(y(:,:,j))); %TODO hard decision, good?

                    for i=1:P.ChannelLength    
                        %TODO reshape(y_rx(i:i+NumberOfChips-1),SeqLen,NumberOfBits/RX); 
                        y_reshape=reshape(y_rx(i:i+NumberOfChips-1),RATE_BITS/6*RX,P.HamLen);
                        [~,indx]=ismember(y_reshape,HadamardMatrix,'rows');
                        rxbits = reshape(de2bi(max(0,indx-1), 6), 1, 516); %TODO magick number + max(0, ..) WTF????
                        
                        % TODO remove this???
                        %rxsymbols(i,:)=SpreadSequence(:,j).'*y_reshape;
                        %rxsymbols(i,:)=rxsymbols(i,:)*conj(himp(j,i));
                    end
                    %TODO: reshape(sum(rxsymbols,1) < 0,1,P.NumberOfSymbols);
                    %rxbits(j:RX:RX*RATE_BITS) = reshape(sum(rxsymbols,1) < 0,1,RATE_BITS);
                end
            otherwise,
                disp('Source Encoding not supported')
        end
        
        % TODO: Decoding Viterbi, is that correct??? Dunno if we have to use it
        % :S
        % rxbits = step(convDec,rxbits.'); right now it doesnt work,
        % problem of shape and size of array
        rxbits = vitdec(rxbits,trellis,34,'trunc','hard');

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
