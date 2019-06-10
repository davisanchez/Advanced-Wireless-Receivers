% Wireless Receivers Project:
% Anael Buchegger, Tim Tuuva, David Sanchez
%
% Direct Sequence Spread Spectrum Simulation Framework
%
% Telecommunications Circuits Laboratory
% EPFL


function BER = MIMOsimulatorII(P)

    RX = P.RXperUser;
    TX = P.TXperUser;
    NUsers = P.CDMAUsers;
    
    % Generate the spreading sequences
    HadamardMatrix = hadamard(P.HadLen)/sqrt(P.HadLen);    %TODO normalization      
    
    NbTXBits    = P.Rate*(P.NumberOfBits + P.Q_Ind + P.K-1);
    NumberOfChips  = NbTXBits*P.HadLen;
    
    HadamSequence     = HadamardMatrix(:,42);
    
    LongCode = comm.PNSequence('Polynomial',[42 7 6 5 3 2 1 0], ...
                               'Mask', P.SequenceMask, ...
                               'InitialConditions', randi([0 1],1,42), ...
                               'SamplesPerFrame', NbTXBits);
                           
    PNSequence = step(LongCode).';
    
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

figure;
for frame = 1:P.NumberOfFrames
    
    frame
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bits = randi([0 1],TX,P.NumberOfBits); % Random Data, High rate
    
    % TODO debug MIMO USE THIS TO DEBUG WITH RX = 2 !!
    % bits = randi([0 1],1,P.NumberOfBits); % Random Data
    % bits = [bits; bits];
    
    % Add Frame Quality Indicator (bonus)
    bits_ind = [bits randi([0 1],TX,P.Q_Ind)];
    
    % Convolutional encoding
    encoded_bits = zeros(TX,NbTXBits);
    for t=1:TX
        encoded_bits(t,:) = convEnc(bits_ind(t,:).').';
    end
    
    % Symbol repetition
    % encoded_bits = repmat(encoded_bits, 1, NbTXBits/length(encoded_bits));
    
    % Interleaver
    if strcmp(P.Interleaving, 'On')
        encoded_bits=matintrlv(encoded_bits.',32,12).';    
    end
    
    % Pulse Shape (PNSequence)
    PN_symbols = zeros(TX, NbTXBits);
    for t=1:TX
        PN_symbols(t,:) = xor(PNSequence, encoded_bits(t,:)); % TODO add user loop here
    end

    % Modulation : BPSK
    symbols = -(2*PN_symbols - 1);
    
    % Spreading with Hadamard
    symbol_spread = zeros(TX,P.HadLen, NbTXBits);
    for t = 1:TX
        symbol_spread(t,:,:) = HadamSequence * symbols(t,:); % TODO MIMO, what shall we do? Different Hadamard sequence?
        % Would be different for each user
    end
    
    waveform = reshape(symbol_spread,TX,NumberOfChips);   % Only for TX right ? 
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Channel
    switch P.ChannelType
        case {'ByPass','AWGN'}
            himp = 1;
        case 'Multipath',
            himp = sqrt(1/2)* (randn(RX,TX,P.ChannelLength) +...
                            1i * randn(RX,TX,P.ChannelLength));
                        
            for r = 1:RX % Normalization for each combination
                for t = 1:TX
                    himp(r,t,:) = himp(r,t,:)/sqrt(sum(abs(himp(r,t,:)).^2)); 
                end
            end
        case 'Fading',
            himp = channel(P.ChannelLength,NumberOfChipsRX,1,P.CoherenceTime,1); %TODO MIMO add a TX channel?           
        otherwise,
            disp('Channel not supported')
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Simulation
    
    switch P.ChannelType
        case {'AWGN','ByPass'}
            snoise = (randn(1,NumberOfChips) + ...
                      1i* randn(1,NumberOfChips) );
        case {'Multipath','Fading'}
            % Duplicate <-> Diversity?? Along the Channels, should it be along
            % Antennas?? TODO TODO IMPORTANT
            snoise = (randn(RX,NumberOfChipsRX) + ...
                      1i* randn(RX,NumberOfChipsRX) );    
            %mwaveform = repmat(waveform,[1 P.ChannelLength 1]);
    end
    
    % SNR Range
    for ss = 1:length(P.SNRRange)
        SNRdb  = P.SNRRange(ss);
        SNRlin = 10^(SNRdb/10);
        noise  = 1/sqrt(2*P.HadLen*SNRlin) *snoise;
        
        % Channel
        switch P.ChannelType
            case 'ByPass', 
                y = waveform * himp;
            case 'AWGN',
                y = waveform * himp + noise;  
            case 'Multipath'
                y = zeros(RX,NumberOfChipsRX); %Normally add the users here!
                for r = 1:RX
                    for t = 1:TX
                        y(r,:) = y(r,:) + conv(squeeze(waveform(t,:)),squeeze(himp(r,t,:)));  
                    end
                end
                y = y + noise;
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
                rxsymbols = conj(himp).*HadamSequence.'*reshape(y, P.HadLen, NbTXBits);
                
                desp_bits = rxsymbols < 0;
                
            case {'Multipath','Fading'} % RAKE receiver (will we have other types ?)
                % Despreading
                rxsymbols = zeros(RX,P.RakeFingers,NbTXBits); % TODO diversity or High rate possibility???
                desp_bits = zeros(TX,NbTXBits);
                
                % Separation between antennas ? How to get himp ?
                for r = 1:RX
                    for finger = 1:P.RakeFingers
                        
                        if strcmp(P.ChannelType,'Multipath')
                            % Despreading
                            rx_despread = HadamSequence.'*reshape(...
                                y(r,finger:finger+NumberOfChips-1),...
                                P.HadLen, NbTXBits);
                        end
                        % Neutralizing global channel effect
                        rxsymbols(r,finger,:) = rx_despread;
                        
                    end
                end
            otherwise
                disp('Source Encoding not supported')
        end
                
        %% MIMO detector
        switch P.Detector
            case 'ZF'
                % Zero Forcing Detector
                % We permute the himp to combine such that the H-1 will
                % cancel out the interference inside rxsymbols
                H = reshape(permute(himp, [1,3,2]), P.ChannelLength*P.RXperUser, P.TXperUser); % TODO need to select the best himp within the rake!
                desp_bits = pinv(H) * reshape(rxsymbols, P.ChannelLength*P.RXperUser, []);
                %plot(desp_bits.', '.')
                
            case 'MMSE'
                % MMSE
                H = reshape(permute(himp, [1,3,2]), P.ChannelLength*P.RXperUser, P.TXperUser); % TODO need to select the best himp within the rake!
                Ps = 1 / (SNRlin);
                G = (H' * H + (P.TXperUser / Ps) * eye(P.TXperUser)) \ H'; %N_T lower rate streams = TXperuser
                desp_bits = G * reshape(rxsymbols, P.ChannelLength*P.RXperUser, []);
                %plot(desp_bits.', '.')
                
                
            case 'SIC'
                H = reshape(permute(himp, [1,2,3]), P.ChannelLength*P.RXperUser, P.TXperUser); % TODO need to select the best himp within the rake!
                H_k = H;
                y_k = reshape(rxsymbols, P.ChannelLength*P.RXperUser, []);
                for k = 1:P.TXperUser
                    
                    e = [1, zeros(1, P.TXperUser - k)];
                    g = e * pinv(H_k);
                    
                    % Estimate signal
                    s_est = (g*y_k);
                    
                    y_k = y_k - H(:, k)*s_est;
                    %plot(y_k.', '.')
                    
                    % Remove k-th column of H_k
                    H_k = H_k(:, 2:end);
                end
                desp_bits = y_k;
            otherwise
                disp('Detector not supported')
        end
        %% Hard Decision
        desc_bits = desp_bits < 0;
        
        %% UN-PN
        unPN_symbols = zeros(TX, NbTXBits);
        for t=1:TX
            unPN_symbols(t,:) = xor(PNSequence, desc_bits(t,:)); % TODO add user loop here
        end
        
        % De-interleaver
        unPN_symbols=double(unPN_symbols);   
        if strcmp(P.Interleaving, 'On')
           unPN_symbols=matdeintrlv(unPN_symbols.',32,12).';
        end    
        
        %% Decoding Viterbi
        decoded_bits = zeros(TX,NbTXBits/2);
        for t = 1:TX
            decoded_bits(t,:) = convDec(unPN_symbols(t,:).').'; % TODO, beurk beurk no??
        end
        
        %% Remove the 8 bit encoding trail
        rxbits = decoded_bits(:,1:end-P.Q_Ind-8); % TODO magick number

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% BER count
        errors =  sum(sum(rxbits ~= bits)); % TODO good way to compute error here?
        Results(ss) = Results(ss) + errors;
        
    end
end

BER = Results/(P.NumberOfBits*P.NumberOfFrames*RX); %TODO added RX here seems logic
end