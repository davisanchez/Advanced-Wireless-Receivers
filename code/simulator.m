% Wireless Receivers Project:
% Anael Buchegger, Tim Tuuva, David Sanchez
%
% Direct Sequence Spread Spectrum Simulation Framework
%
% Telecommunications Circuits Laboratory
% EPFL


function BER = simulator(P)

    RX = P.RXperUser;
    TX = P.TXperUser;
    NUsers = P.CDMAUsers;
    
    % Generate the spreading sequences (One per antenna ?)
    HadamardMatrix = hadamard(P.HadLen)/sqrt(P.HadLen);    %TODO normalization      
    
    % 
    NumberOfBits   = P.NumberOfBits*P.Modulation; % per Frame
    
    %RATE_BITS = P.Rate * NumberOfBits;
    NbTXBits    = P.Rate*(P.NumberOfBits + P.Q_Ind + P.K-1);
    NumberOfChips  = NbTXBits*P.HadLen;
    
    HadamSequence     = HadamardMatrix(:,42);
    
    LongCode = comm.PNSequence('Polynomial',[42 7 6 5 3 2 1 0], ...
                               'Mask', P.SequenceMask, ...
                               'InitialConditions', randi([0 1],1,42), ...
                               'SamplesPerFrame', NbTXBits);
    % TODO MIMO ??? is that ok?
    PNSequence = zeros(NUsers,NbTXBits); % SHOULD BE CDMA USERS related no??
    for n = 1:NUsers %each user has a PN
        PNSequence(n,:) = step(LongCode);
    end
    
    % Channel
    switch P.ChannelType
        case 'Multipath'
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
    bits = randi([0 1],TX,P.NumberOfBits); % Random Data
    
    % TODO debug MIMO USE THIS TO DEBUG WITH RX = 2 !!
    %bits = randi([0 1],1,NumberOfBits); % Random Data
    %bits = [bits; bits];
    
    % Add Frame Quality Indicator (bonus)
    bits_ind = [bits randi([0 1],TX,P.Q_Ind)];
    
    % Convolutional encoding
    encoded_bits = zeros(TX,NbTXBits);
    for t=1:TX
        encoded_bits(t,:) = convEnc(bits_ind(t,:).').'; % TODO doesnt give the same encoding signal, why??? IMPORTANT
    end
    
    % Symbol repetition
    encoded_bits = repmat(encoded_bits, 1, NbTXBits/length(encoded_bits));
    
    % Here comes the interleaver (TODO) Also need de-interleaver !
%     encoded_bits=encoded_bits.';
%     encoded_bits(:,1)=matintrlv(encoded_bits(:,1),32,12);
%     encoded_bits(:,2)=matintrlv(encoded_bits(:,2),32,12);
%     encoded_bits=encoded_bits.';
    
    % Pulse Shape (PNSequence)
    
    PN_symbols = zeros(TX, NbTXBits);
    %user loop
    for k=1:NUsers
        for t=1:TX
            PN_symbols(t,:) = xor(PNSequence(k,:), encoded_bits(t,:)); % TODO add user loop here
        end
    end

    % Modulation : BPSK
    symbols = -(2*PN_symbols - 1);
    
    % Spreading with Hadamard
    symbol_spread = zeros(TX,P.HadLen, NbTXBits);
    for t = 1:TX
        symbol_spread(t,:,:) = HadamSequence * symbols(t,:); % TODO MIMO, what shall we do? Different Hadamard sequence?
    end
    
    waveform = reshape(symbol_spread,TX,1,P.HadLen*NbTXBits);   % Only for TX right ? 
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Channel
    switch P.ChannelType
        case {'ByPass','AWGN'}
            himp = ones(RX,TX,P.ChannelLength,1);
        case 'Multipath',
            himp = sqrt(1/2)* (randn(RX,TX,P.ChannelLength,1) + 1i * randn(RX,TX,P.ChannelLength,1));
        case 'Fading',
            himp = channel(P.ChannelLength,NumberOfChipsRX,1,P.CoherenceTime,1); %TODO MIMO add a TX channel?
        
        otherwise,
            disp('Channel not supported')
    end
    
    % Duplicate <-> Diversity?? Along the Channels, should it be along
    % Antennas?? TODO TODO IMPORTANT
    mwaveform = repmat(waveform,[1 P.ChannelLength 1]);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Simulation
    snoise = ( randn(RX,TX,P.ChannelLength,NumberOfChips) + 1i* randn(RX,TX,P.ChannelLength,NumberOfChips) );
    
    % SNR Range
    for ss = 1:length(P.SNRRange)
        SNRdb  = P.SNRRange(ss);
        SNRlin = 10^(SNRdb/10);
        noise  = 1/sqrt(2*P.HadLen*SNRlin) *snoise;
        
        % Channel
        % y -> (Users, Antenna, h_Channel, Data) 
        y = zeros(RX,P.ChannelLength,NumberOfChipsRX); %Normally add the users here!
        switch P.ChannelType
            case 'ByPass',
                for r = 1:RX
                    for t = 1:TX
                        for i = 1:P.ChannelLength
                            % This is needed because Matlab cannot handle conv
                            % with matrix like structure even though it is a
                            % vector
                            % We sum all the contribution from every TX
                            % antennas !!! and Also sum from every user 
                            % In real life we can't separate the data as a
                            % row of a given matrix XP tuff life :-/
                            signal = squeeze(mwaveform(t,i,:));
                            y(r,i,:) = squeeze(y(r,i,:)) + ...
                                       squeeze(conv(signal,himp(r,t,i,:)));
                        end
                    end
                end
            case 'AWGN',
                for r = 1:RX
                    for t = 1:TX
                        for i = 1:P.ChannelLength
                            signal = squeeze(mwaveform(t,i,:));
                            y(r,i,:) = squeeze(y(r,i,:)) + ...
                                       squeeze(conv(signal,himp(r,t,i,:))) + ...
                                       squeeze(noise(r,t,i,:));
                        end
                    end
                end
            case 'Multipath'
                for r = 1:RX
                    for t = 1:TX
                        for i = 1:P.ChannelLength
                            signal = squeeze(mwaveform(r,i,:));
                            y(r,i,i:i+NumberOfChips-1) = squeeze(y(r,i,i:i+NumberOfChips-1)) + ...
                                                         squeeze(conv(signal,himp(r,t,i,:))) + ...
                                                         squeeze(noise(r,t,i,:));
                        end
                    end
                end
            case 'Fading',
                for i = 1:P.ChannelLength
                    y(i,:,RrX) = mwaveform(i,:,r) .* himp(i,:,r) + noise;
                end
            otherwise,
                disp('Channel not supported')
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Receiver
        switch P.ReceiverType
            case 'Rake',  
                % Despreading
                rxsymbols = zeros(RX,P.RakeFingers,P.ChannelLength*NbTXBits); % TODO diversity or High rate possibility???

                
                for r = 1:RX
                    if ~strcmp(P.ChannelType,'Fading')
                        [~,ind] = maxk(himp(r,:,:,:),P.ChannelLength,3);
                    else
                        [himp_mean,ind] = maxk(mean(himp(r,:,:)),P.ChannelLength);
                    end
                    for finger = 1:P.RakeFingers
                        if ~strcmp(P.ChannelType,'Fading')
                            
                            rx_channel = squeeze(conj(himp(r,ind(finger),:)));
                            rx_despread = HadamSequence.'*reshape(...
                                                          y(r,ind(finger),ind(finger):ind(finger)+NumberOfChips-1),...
                                                          P.HadLen, NumberOfChips/P.HadLen);

                            rxsymbols(r,finger,:) = reshape(rx_channel * rx_despread, 1, P.ChannelLength*NbTXBits);
                        else
                            rxsymbols(r,finger,:) = squeeze(conj(himp_mean(r,ind(finger),:)))*HadamSequence.'*...
                                reshape(y(r,ind(finger),:), ...
                                P.HadLen, NumberOfChipsRX/P.HadLen);
                        end
                    end
                end
                % TODO diversity or High rate possibility???
                desc_bits = zeros(RX,NbTXBits);
                for r = 1:RX
                    % THIS IS AN EXAMPLE OF DIVERSITY ( I guess... ) TODO
                    % TODO IMPORTANT HERE, the sum is doing the diversity
                    desc_bits(r,:) = sum(reshape(sum(rxsymbols(r,:,:),2) < 0,P.ChannelLength,NbTXBits),1)/P.ChannelLength; 
                end
                
            otherwise,
                disp('Source Encoding not supported')
        end
        % UN-PN
        unPN_symbols = zeros(RX, NbTXBits);
        for r=1:RX
            unPN_symbols(r,:) = xor(PNSequence(NUsers,:), desc_bits(r,:)); % TODO add user loop here
        end
        
        
        % Decoding Viterbi
        decoded_bits = zeros(RX,NbTXBits/2);
        for r = 1:RX
            decoded_bits(r,:) = convDec(double(unPN_symbols(r,:)).').'; % TODO, beurk beurk no??
        end
        
        % Remove the 8 bit encoding trail
        rxbits = decoded_bits(:,1:end-P.Q_Ind-8); % TODO magick number

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % BER count
        errors =  sum(sum(rxbits ~= bits)); % TODO good way to compute error here?
        
        Results(ss) = Results(ss) + errors;
        
    end
end

BER = Results/(P.NumberOfBits*P.NumberOfFrames*RX); %TODO added RX here seems logic
end