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
                           
    PNSequence = step(LongCode);
    
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
    bits = randi([0 1],P.NumberOfBits); % Random Data
    
    % Add Frame Quality Indicator (bonus)
    bits_ind = [bits randi([0 1],TX,P.Q_Ind)];
    
    % Convolutional encoding
    encoded_bits = convEnc(bits_ind.').'; % TODO doesnt give the same encoding signal, why??? IMPORTANT
    
    % Symbol repetition
    encoded_bits = repmat(encoded_bits, 1, NbTXBits/length(encoded_bits));
    
    % Here comes the interleaver (TODO)
    
    % Pulse Shape (PNSequence)
    PN_symbols = zeros(TX, NbTXBits);
    for t=1:TX
        PN_symbols(t,:) = xor(PNSequence(NUsers,:), encoded_bits(t,:)); % TODO add user loop here
    end

    % Modulation : BPSK
    symbols = -(2*PN_symbols - 1);
    
    % Spreading with Hadamard
    symbol_spread = zeros(TX,P.HadLen, NbTXBits);
    for t = 1:TX
        symbol_spread(t,:,:) = HadamSequence * symbols(t,:); % TODO MIMO, what shall we do? Different Hadamard sequence?
    end
    
    waveform = reshape(symbol_spread,TX,P.HadLen*NbTXBits);   % Only for TX right ? 
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Channel
    switch P.ChannelType
        case {'ByPass','AWGN'}
            himp = ones(RX,TX);
        case 'Multipath',
            himp = sqrt(1/2)* (randn(RX,TX,P.ChannelLength,1) + 1i * randn(RX,TX,P.ChannelLength,1));
        case 'Fading',
            himp = channel(P.ChannelLength,NumberOfChipsRX,1,P.CoherenceTime,1); %TODO MIMO add a TX channel?           
        otherwise,
            disp('Channel not supported')
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Simulation
    
    switch P.ChannelType
        case {'AWGN','ByPass'}
            snoise = (randn(RX,TX,NumberOfChips) + ...
                      1i* randn(RX,TX,NumberOfChips) );
        case {'Multipath','Fading'}
            % Duplicate <-> Diversity?? Along the Channels, should it be along
            % Antennas?? TODO TODO IMPORTANT
            snoise = (randn(RX,TX,P.ChannelLength,NumberOfChips) + ...
                      1i* randn(RX,TX,P.ChannelLength,NumberOfChips) );    
            mwaveform = repmat(waveform,[1 P.ChannelLength 1]);
    end
    
    % SNR Range
    for ss = 1:length(P.SNRRange)
        SNRdb  = P.SNRRange(ss);
        SNRlin = 10^(SNRdb/10);
        noise  = 1/sqrt(2*P.HadLen*SNRlin) *snoise;
        
        % Channel
        switch P.ChannelType
            case 'ByPass',
                % y -> (Users, Antenna, h_Channel, Data) 
                y = zeros(RX,NumberOfChipsRX); %Normally add the users here!
                for r = 1:RX
                    for t = 1:TX
                        % This is needed because Matlab cannot handle conv
                        % with matrix like structure even though it is a
                        % vector
                        % We sum all the contribution from every TX
                        % antennas !!! and Also sum from every user 
                        % In real life we can't separate the data as a
                        % row of a given matrix XP tuff life :-/
                        signal = squeeze(waveform(t,:));
                        y(r,:) = squeeze(y(r,:)) + ...
                                   squeeze(conv(signal,himp(r,t,:)));
                    end
                end
            case 'AWGN',
                % y -> (Users, Antenna, h_Channel, Data) 
                y = zeros(RX,NumberOfChipsRX); %Normally add the users here!
                for r = 1:RX
                    for t = 1:TX
                        signal = waveform(t,:);
                        y(r,:) = y(r,:)+ signal*himp(r,t) + squeeze(noise(r,t,:)).';
                    end
                end
            case 'Multipath'
                % y -> (Users, Antenna, h_Channel, Data) 
                y = zeros(RX,P.ChannelLength,NumberOfChipsRX); %Normally add the users here!
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
                % y -> (Users, Antenna, h_Channel, Data) 
                y = zeros(RX,P.ChannelLength,NumberOfChipsRX); %Normally add the users here!
                for i = 1:P.ChannelLength
                    y(i,:,RrX) = mwaveform(i,:,r) .* himp(i,:,r) + noise;
                end
            otherwise,
                disp('Channel not supported')
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Receiver
        switch P.ChannelType % Rake receiver makes no sense in this case
            case {'AWGN','ByPass'}
                % Despreading
                rxsymbols = zeros(RX,NbTXBits); % TODO diversity or High rate possibility???
                desc_bits = zeros(RX,NbTXBits);
                    
                for r = 1:RX
                    rx_channel = conj(himp(r,:)); % Not quite sure what to do with himp here
                    rxsymbols(r,:) = HadamSequence.'*reshape(y(r,:), P.HadLen, NbTXBits);
                    
                    % THIS IS AN EXAMPLE OF DIVERSITY ( I guess... ) TODO
                    % TODO IMPORTANT HERE, the sum is doing the diversity
                    desc_bits(r,:) = sum(reshape(sum(rxsymbols(r,:),2) < 0,NbTXBits),1); 
                end
            
            case {'Multipath','Fading'} % RAKE receiver (will we have other types ?)
                % Despreading
                rxsymbols = zeros(RX,P.RakeFingers,P.ChannelLength*NbTXBits); % TODO diversity or High rate possibility???
                
                for r = 1:RX
                    % Order the best fingers
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

BER = Results/(NumberOfBits*P.NumberOfFrames*RX); %TODO added RX here seems logic
end